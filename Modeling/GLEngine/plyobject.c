#include "plyobject.h"
#include "utils.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/glext.h>
#else
#include <GL/glut.h>
#include <GL/glext.h>
#endif

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>

PLYObject *readPLYFile(char *filename)
{
  PLYObject *plyObject = malloc(sizeof(PLYObject));

  FILE *inputFile = fopen(filename,"r");

  char *lineBuffer = malloc(sizeof(char) * MAXLINELENGTH);
  char *tokenBuffer = malloc(sizeof(char) * MAXTOKENLENGTH);
  
  bool readingHeader = true;
  bool readingVertices = false;
  bool readingFaces = false;
  int numVerticesRead = 0;
  int numFacesRead = 0;
  int i;

  printf("Attempting to Load mesh from .ply file: %s\n",filename);
  
  while(fgets(lineBuffer,MAXLINELENGTH,inputFile)!= NULL)
    {
      tokenBuffer = strtok(lineBuffer," \t\n");
      if (strcmp(tokenBuffer,"end_header") == 0)
        {
          readingHeader = false;
          readingVertices = true;

          //Allocate space for vertices
          plyObject->_vertices = malloc(sizeof(GLfloat) * plyObject->_numVertices * 3);
        }
      else if (strcmp(tokenBuffer,"element") == 0)
        {
          tokenBuffer = strtok(NULL," \t\n");
          if (strcmp(tokenBuffer,"vertex") == 0)
            {
              tokenBuffer = strtok(NULL," \t\n");
              plyObject->_numVertices = atoi(tokenBuffer);
            }
          else if(strcmp(tokenBuffer,"face") == 0)
            {
              tokenBuffer = strtok(NULL," \t\n");
              plyObject->_numFaces = atoi(tokenBuffer);
            }
        }
      else if(readingHeader == false)
        {
          //Check if we need to start reading faces or end
          if (numVerticesRead == plyObject->_numVertices && readingVertices == true)
            {
              readingVertices = false;
              readingFaces = true;

              //allocate space for faces
              plyObject->_numVertPerFace = atoi(tokenBuffer);
              plyObject->_faces = malloc(sizeof(int) * plyObject->_numVertPerFace * plyObject->_numFaces);
            }

          if (numFacesRead == plyObject->_numFaces)
            {
              fprintf(stderr,"Error: About to read one more face that was stated...\n");
            }
          
          if (readingVertices == true)
            {
              for (i=0; i<3; i++)
                {
                  plyObject->_vertices[numVerticesRead * 3 + i] = atof(tokenBuffer);
                  tokenBuffer = strtok(NULL," \t\n");
                }
              numVerticesRead++;
            }
          else if (readingFaces == true)
            {
              if (atoi(tokenBuffer) != plyObject->_numVertPerFace)
                {
                  fprintf(stderr,"Error: Faces did not have consistant number of verticies\n");
                  exit(1);
                }
              for (i=0; i<plyObject->_numVertPerFace; i++)
                {
                  tokenBuffer = strtok(NULL," \t\n");
                  plyObject->_faces[numFacesRead * plyObject->_numVertPerFace + i] = atof(tokenBuffer);
                }
              numFacesRead++;
            }
        }
    }

  //Check if you properly read in a ply file otherwise reject the plyfile
  if (numFacesRead != plyObject->_numFaces || numVerticesRead != plyObject->_numVertices)
    {
      fprintf(stderr,"Error: Number of vertices and faces to read in do not match file\n");
      exit(1);
    }

  printf("Mesh loaded successfully\n");
  
  CalcPLYObjectNormals(plyObject,true);
  
  return plyObject;
}

void NormalizePLYObject(PLYObject *p)
{
  assert(p->_numVertices > 0);
  
  int i;

  float minX,minY,minZ,maxX,maxY,maxZ;

  minX = p->_vertices[0]; maxX = p->_vertices[0];
  minX = p->_vertices[1]; maxX = p->_vertices[1];
  minX = p->_vertices[2]; maxX = p->_vertices[2];

  //Find the min/max for x,y,z
  for (i=0; i<p->_numVertices; i++)
    {
      if (p->_vertices[i*3] < minX)
        {
          minX = p->_vertices[i*3];
        }
      if (p->_vertices[i*3] > maxX)
        {
          maxX = p->_vertices[i*3];
        }
      if (p->_vertices[i*3+1] < minY)
        {
          minY = p->_vertices[i*3+1];
        }
      if (p->_vertices[i*3+1] > maxY)
        {
          maxY = p->_vertices[i*3+1];
        }
      if (p->_vertices[i*3+2] < minZ)
        {
          minZ = p->_vertices[i*3+2];
        }
      if (p->_vertices[i*3+2] > maxZ)
        {
          maxZ = p->_vertices[i*3+2];
        }
    }

  float scaleFactor = 1 / (max(maxX - minZ, max(maxY - minY, maxZ - minZ)));
  //Uniform scaling such that each point lies within 0 - 1
  for (i=0; i<p->_numVertices; i++)
    {
      p->_vertices[3*i] = scaleFactor*(p->_vertices[3*i] - minX);
      p->_vertices[3*i+1] = scaleFactor*(p->_vertices[3*i+1] - minY);
      p->_vertices[3*i+2] = scaleFactor*(p->_vertices[3*i+2] - minZ);
    }
}

void CalcPLYObjectNormals(PLYObject *p, bool flipSign)
{
  int i,j;
  GLuint *face;
  GLfloat *norm;
  GLfloat len;
  GLfloat recip;
  float x,y,z;
  float x0,y0,z0;
  float x1,y1,z1;

  p->_numNormals = p->_numVertices;
  p->_normals = malloc(sizeof(GLfloat) * p->_numNormals * 3);
  
  //Zero out all the normals
  for (i=0; i<(p->_numNormals * 3); i++)
    {
      p->_normals[i] = 0;
    }

  // Find the normal of each face , add it to each vertex adjacent
  // To the face.
  for (i=0; i< p->_numFaces; i++)
    {
      face = &(p->_faces[i*p->_numVertPerFace]);

      // Finding two vectors that represent the plane
      x0 = p->_vertices[face[2]*3] - p->_vertices[face[0]*3];
      y0 = p->_vertices[face[2]*3+1] - p->_vertices[face[0]*3+1];
      z0 = p->_vertices[face[2]*3+2] - p->_vertices[face[0]*3+2];

      x1 = p->_vertices[face[1]*3] - p->_vertices[face[0]*3];
      y1 = p->_vertices[face[1]*3+1] - p->_vertices[face[0]*3+1];
      z1 = p->_vertices[face[1]*3+2] - p->_vertices[face[0]*3+2];

      /* find cross-product between these vectors */
      x = y0 * z1 - z0 * y1;
      y = z0 * x1 - x0 * z1;
      z = x0 * y1 - y0 * x1;

      /* normalize this vector */
      len = x*x + y*y + z*z;
      if (len == 0) {
        x = y = z = 0;
      }
      else {
        recip = 1 / sqrt (len);
        x *= recip;
        y *= recip;
        z *= recip;
      }

      /* add this normal to each vertex that is adjacent to face */
      for (j = 0; j < p->_numVertPerFace; j++) {
        p->_normals[face[j]*3] += x;
        p->_normals[face[j]*3+1] += y;
        p->_normals[face[j]*3+2] += z;
      }
    }

   /* normalize all the normals at the vertices */

  for (i = 0; i < p->_numNormals; i++) {
    norm = &(p->_normals[i*3]);
    len = norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2];
    if (len == 0) {
      norm[0] = 0;
      norm[1] = 0;
      norm[2] = 0;
    }
    else {
      if (flipSign)
        recip = -1 / sqrt (len);
      else
        recip = 1 / sqrt (len);
      norm[0] *= recip;
      norm[1] *= recip;
      norm[2] *= recip;
    }
  }
}

/*
 *  Sends PLYObject vertices, normals, and faces
 *  to the GPU. Returns integer values specifying
 *  where the buffers for each one is located
 *  this allows for multiple objects to be created
 */
GLuint *initializePLYObject(PLYObject *p)
{
  GLuint *bufObjects = malloc(sizeof(GLuint)*3);

  glGenBuffers(3,bufObjects);
  // 0 stores vertices
  // 1 stores normals
  // 2 stores indices
  // These are used to grab the buffer data that we need

  //Send vertex coordinates to the GPU
  glBindBuffer(GL_ARRAY_BUFFER,bufObjects[0]);
  glBufferData(GL_ARRAY_BUFFER,sizeof(GLfloat)*3*p->_numVertices,p->_vertices,GL_STATIC_DRAW);
  glVertexPointer(3,GL_FLOAT,0,NULL);
  glEnableClientState(GL_VERTEX_ARRAY); 
  
  //Send normal coordinates to the GPU
  glBindBuffer(GL_ARRAY_BUFFER,bufObjects[1]);
  glBufferData(GL_ARRAY_BUFFER,sizeof(GLfloat)*3*p->_numNormals,p->_normals,GL_STATIC_DRAW);
  glNormalPointer(GL_FLOAT,0,NULL);
  glEnableClientState(GL_NORMAL_ARRAY);

  //Setup the faces
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,bufObjects[2]);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(GLuint)*p->_numVertPerFace*p->_numFaces,p->_faces,GL_STATIC_DRAW);
  //glIndexPointer(GL_UNSIGNED_INT,0,NULL);
  //glEnableClientState(GL_INDEX_ARRAY);
  
  checkGLErrors("Loading ply mesh");
  return bufObjects;
}

void DrawPLYObject(PLYObject *p)
{
  

}


/*
 *  Prints the data contained within the PlyObject
 *  Use for relativly small ply files otherwise
 *  this will be HUGE. Make sure to call after normals
 *  have already been formed or it will SEGFAULT
 */
void PrintPLYObject(PLYObject *p)
{
  int i,j;
  printf("The Vertices: %d\n",p->_numVertices);
  for (i=0; i<p->_numVertices; i++)
    {
      printf("%f %f %f\n",p->_vertices[3*i],p->_vertices[1+3*i],p->_vertices[2+3*i]);
    }
  
  printf("The Normals: %d\n",p->_numNormals);
  for (i=0; i<p->_numNormals; i++)
    {
      printf("%f %f %f\n",p->_normals[3*i],p->_normals[1+3*i],p->_normals[2+3*i]);
    }
  
  printf("The Faces: %d\n",p->_numFaces);
  for (i=0; i<p->_numFaces; i++)
    {
      for (j=0; j<p->_numVertPerFace; j++)
        {
          printf("%d ",p->_faces[j+p->_numVertPerFace*i]);
        }
      printf("\n");
    }
}
