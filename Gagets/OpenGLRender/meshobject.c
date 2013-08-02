#include "meshobject.h"
#include "utils.h"
#include "texture.h"
#include "lib/glmath.h"
#include "shader.h"

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#include <OpenGL/glext.h>
#else
#include <GL/gl.h>
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

MeshObject *readPLYFile(char *filename)
{
  MeshObject *meshObject = malloc(sizeof(MeshObject));

  FILE *inputFile = fopen(filename,"r");

  if (inputFile == NULL)
    {
      fprintf(stderr, "Error: Unable to Load Mesh from .ply file: %s\n",filename);
      free(meshObject);
      return NULL;
    }
  
  char *lineBuffer = malloc(sizeof(char) * MAXLINELENGTH);
  char *tokenBuffer = malloc(sizeof(char) * MAXTOKENLENGTH);
  
  bool readingHeader = true;
  bool readingVertices = false;
  bool readingFaces = false;
  int numVerticesRead = 0;
  int numFacesRead = 0;
  int maxAllocationSize;
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
          meshObject->_vertices = malloc(sizeof(vertex) * meshObject->_numVertices);
        }
      else if (strcmp(tokenBuffer,"element") == 0)
        {
          tokenBuffer = strtok(NULL," \t\n");
          if (strcmp(tokenBuffer,"vertex") == 0)
            {
              tokenBuffer = strtok(NULL," \t\n");
              meshObject->_numVertices = atoi(tokenBuffer);
            }
          else if(strcmp(tokenBuffer,"face") == 0)
            {
              tokenBuffer = strtok(NULL," \t\n");
              meshObject->_numFaces = atoi(tokenBuffer);
            }
        }
      else if(readingHeader == false)
        {
          //Check if we need to start reading faces or end
          if (numVerticesRead == meshObject->_numVertices && readingVertices == true)
            {
              readingVertices = false;
              readingFaces = true;

              //allocate space for faces
              meshObject->_faces = malloc(sizeof(face) * meshObject->_numFaces);
              maxAllocationSize = meshObject->_numFaces;
            }

          if (numFacesRead == meshObject->_numFaces)
            {
              fprintf(stderr,"Error: About to read one more face than was stated...\n");
            }
          
          if (readingVertices == true)
            {
              for (i=0; i<3; i++)
                {
                  meshObject->_vertices[numVerticesRead]._xyzCoord[i] = atof(tokenBuffer);
                  tokenBuffer = strtok(NULL," \t\n");
                }
              numVerticesRead++;
            }
          else if (readingFaces == true)
            {
              int numVerticesForFace = atoi(tokenBuffer);
              int *indices = malloc(sizeof(int) * numVerticesForFace);
              for (i=0; i<numVerticesForFace; i++)
                {
                  tokenBuffer = strtok(NULL," \t\n");
                  indices[i] = atoi(tokenBuffer);
                }

              // Might Have to grow the size of the array
              // Because we have to convert all faces to triangles
              if (numFacesRead + numVerticesForFace - 2 > maxAllocationSize)
                {
                  meshObject->_faces = realloc(meshObject->_faces,sizeof(face)*maxAllocationSize*2);
                  maxAllocationSize *= 2;
                  if (meshObject->_faces == NULL)
                    {
                      fprintf(stderr,"Error: failed to reallocate memory for faces\n");
                    }
                }
              
              //Assumes that polygons read in are convex
              for (i=0; i<numVerticesForFace-2; i++)
                {
                  meshObject->_faces[numFacesRead + i]._indices[0] = indices[0];
                  meshObject->_faces[numFacesRead + i]._indices[1] = indices[i+1];
                  meshObject->_faces[numFacesRead + i]._indices[2] = indices[i+2];
                }

              free(indices);
              
              meshObject->_numFaces += numVerticesForFace - 3;
              numFacesRead += numVerticesForFace - 2;
            }
        }
    }
  meshObject->_faces = realloc(meshObject->_faces, sizeof(face) * meshObject->_numFaces);

  //Check if you properly read in a ply file otherwise reject the plyfile
  if (numFacesRead != meshObject->_numFaces || numVerticesRead != meshObject->_numVertices)
    {
      fprintf(stderr,"Error: Number of vertices and faces to read in do not match file\n");
      exit(1);
    }

  printf("Mesh loaded successfully\n");

  CalcMeshObjectVertexVectors(meshObject, true);
  
  return meshObject;
}

MeshObject *readOBJFile(char *filename)
{
  //Initialize mesh object
  MeshObject *meshObject = malloc(sizeof(MeshObject));
  meshObject->_numFaces = 1;
  meshObject->_faces = malloc(sizeof(face) * 8);

  FILE *inputFile = fopen(filename,"r");

  if (inputFile == NULL)
    {
      fprintf(stderr, "Error: Unable to Load Mesh from .obj file: %s\n",filename);
      free(meshObject);
      return NULL;
    }
  
  char *lineBuffer = malloc(sizeof(char) * MAXLINELENGTH);
  char *tokenBuffer = malloc(sizeof(char) * MAXTOKENLENGTH);

  double *vertices = (double *) malloc(sizeof(double) * 3);
  double *texCoord = (double *) malloc(sizeof(double) * 3);
  double *vertexNorm = (double *) malloc(sizeof(double) * 3);
  double *vertexTan = (double *) malloc(sizeof(double) * 3);
  double *vertexBiTan = (double *) malloc(sizeof(double) * 3);
  
  int i,j;
  int numVerticesRead = 0, numVerticesAlloc = 1;
  int numTexCoordsRead = 0, numTexCoordAlloc = 1;
  int numVertexNormsRead = 0, numVertexNormsAlloc = 1;
  int numVertexTansRead = 0, numVertexTansAlloc = 1;
  int numVertexBiTansRead = 0, numVertexBiTansAlloc = 1;

  int numFacesRead = 0, numFacesAlloc = 8;

  bool readingFaces = false;
  
  printf("Attempting to Load mesh from .obj file: %s\n",filename);
  
  while (fgets(lineBuffer, MAXLINELENGTH, inputFile) != NULL)
    {
      tokenBuffer = strtok(lineBuffer," \n\r\t");

      if (strcmp(tokenBuffer, "mtllib") == 0)
        {
          tokenBuffer = strtok(NULL, " \n\t\r");
          Material *newMaterial = readMaterialFile(tokenBuffer);
          addMaterialToLibrary(newMaterial);
        }
      else if (strcmp(tokenBuffer, "usemtl") == 0)
        {
          tokenBuffer = strtok(NULL, " \n\t\r");
          meshObject->_material = getMaterial(tokenBuffer);
        }
      else if (strcmp(tokenBuffer, "v") == 0)
        {
          if (numVerticesRead == numVerticesAlloc)
            {
              vertices = realloc(vertices, sizeof(double) * 3 * numVerticesAlloc * 2);
              numVerticesAlloc *= 2;
            }

          for (i=0; i<3; i++)
            {
              tokenBuffer = strtok(NULL," \t\n");
              vertices[3*numVerticesRead+i] = atof(tokenBuffer);
            }

          numVerticesRead++;
        }
      else if (strcmp(tokenBuffer, "vt") == 0)
        {
          if (numTexCoordsRead == numTexCoordAlloc)
            {
              texCoord = realloc(texCoord, sizeof(double) * 3 * numTexCoordAlloc * 2);
              numTexCoordAlloc *= 2;
            }

          for (i=0; i<2; i++)
            {
              tokenBuffer = strtok(NULL," \t\n");
              texCoord[2*numTexCoordsRead+i] = atof(tokenBuffer);
            }

          numTexCoordsRead++;
        }
      else if (strcmp(tokenBuffer, "vn") == 0)
        {
          if (numVertexNormsRead == numVertexNormsAlloc)
            {
              vertexNorm = realloc(vertexNorm, sizeof(double) * 3 * numVertexNormsAlloc * 2);
              numVertexNormsAlloc *= 2;
            }

          for (i=0; i<3; i++)
            {
              tokenBuffer = strtok(NULL," \t\n");
              vertexNorm[3*numVertexNormsRead+i] = atof(tokenBuffer);
            }

          numVertexNormsRead++;
        }
      else if (strcmp(tokenBuffer, "vx") == 0)
        {
          if (numVertexTansRead == numVertexTansAlloc)
            {
              vertexTan = realloc(vertexTan, sizeof(double) * 3 * numVertexTansAlloc * 2);
              numVertexTansAlloc *= 2;
            }

          for (i=0; i<3; i++)
            {
              tokenBuffer = strtok(NULL," \t\n");
              vertexTan[3*numVertexTansRead+i] = atof(tokenBuffer);
            }

          numVertexTansRead++;
        }
      else if (strcmp(tokenBuffer, "vy") == 0)
        {
          if (numVertexBiTansRead == numVertexBiTansAlloc)
            {
              vertexBiTan = realloc(vertexBiTan, sizeof(double) * 3 * numVertexBiTansAlloc * 2);
              numVertexBiTansAlloc *= 2;
            }

          for (i=0; i<3; i++)
            {
              tokenBuffer = strtok(NULL," \t\n");
              vertexBiTan[3*numVertexBiTansRead+i] = atof(tokenBuffer);
            }

          numVertexBiTansRead++;
        }
      //It is assumed that all verticies tan norms etc are read in before faces
      //Thus all the vectors v,vt,vn,vx,vy are fully read in. This does most of the
      // dirty work
      else if (strcmp(tokenBuffer, "f") == 0)
        {
          if (readingFaces == false)
            {
              meshObject->_numVertices = numVerticesRead;
              meshObject->_vertices = malloc(sizeof(vertex) * meshObject->_numVertices);
            }

          // Parse the line that stores a face made of indicies
          // structure of each indice is v/vt/vn
          // only v/vt are used since all other parts are calculated
          // by our code. Max number of indicies is 10
          int numIndicesInFace = 0;
          int *indicesInfo = malloc(sizeof(int) * 30);
          while ((tokenBuffer = strtok(NULL," \t\n/")) != NULL)
            {
              indicesInfo[3*numIndicesInFace] = atoi(tokenBuffer);
              tokenBuffer = strtok(NULL," \t\n/");
              indicesInfo[3*numIndicesInFace+1] = atoi(tokenBuffer);
              tokenBuffer = strtok(NULL," \t\n/");
              indicesInfo[3*numIndicesInFace+2] = atoi(tokenBuffer);

              numIndicesInFace++;
            }
          
          //When reading in the faces each vertex will also be created
          if (numFacesRead + numIndicesInFace - 2 > numFacesAlloc)
            {
              meshObject->_faces = realloc(meshObject->_faces, sizeof(face) * numFacesAlloc * 2);
              numFacesAlloc *= 2;
            }

          //Assumes that polygons read in are convex and in order
          //This assumption may be wrong...
          for (i=0; i<numIndicesInFace-2; i++)
            {
              //Set the faces
              meshObject->_faces[numFacesRead + i]._indices[0] = indicesInfo[0]-1;
              meshObject->_faces[numFacesRead + i]._indices[1] = indicesInfo[3*(i+1)]-1;
              meshObject->_faces[numFacesRead + i]._indices[2] = indicesInfo[3*(i+2)]-1;
              //Set the associated vertices info
              for (j=0; j<3; j++)
                {
                  meshObject->_vertices[indicesInfo[0]-1]._xyzCoord[j] = vertices[3*(indicesInfo[0]-1)+j];
                  meshObject->_vertices[indicesInfo[0]-1]._normal[j] = vertexNorm[3*(indicesInfo[2]-1)+j];
                  meshObject->_vertices[indicesInfo[0]-1]._tangent[j] = vertexTan[3*(indicesInfo[0]-1)+j];
                  meshObject->_vertices[indicesInfo[0]-1]._bitangent[j] = vertexBiTan[3*(indicesInfo[0]-1)+j];
                }
              meshObject->_vertices[indicesInfo[0]-1]._uvCoord[0] = texCoord[2*(indicesInfo[1]-1)];
              meshObject->_vertices[indicesInfo[0]-1]._uvCoord[1] = texCoord[2*(indicesInfo[1]-1)+1];
	      
              for (j=0; j<3; j++)
                {
                  meshObject->_vertices[indicesInfo[3*(i+1)]-1]._xyzCoord[j] = vertices[3*(indicesInfo[3*(i+1)]-1)+j];
                  meshObject->_vertices[indicesInfo[3*(i+1)]-1]._normal[j] = vertexNorm[3*(indicesInfo[3*(i+1)+2]-1)+j];
                  meshObject->_vertices[indicesInfo[3*(i+1)]-1]._tangent[j] = vertexTan[3*(indicesInfo[3*(i+1)]-1)+j];
                  meshObject->_vertices[indicesInfo[3*(i+1)]-1]._bitangent[j] = vertexBiTan[3*(indicesInfo[3*(i+1)]-1)+j];
                }
              meshObject->_vertices[indicesInfo[3*(i+1)]-1]._uvCoord[0] = texCoord[2*(indicesInfo[3*(i+1)+1]-1)];
              meshObject->_vertices[indicesInfo[3*(i+1)]-1]._uvCoord[1] = texCoord[2*(indicesInfo[3*(i+1)+1]-1)+1];

	      
              for (j=0; j<3; j++)
                {
                  meshObject->_vertices[indicesInfo[3*(i+2)]-1]._xyzCoord[j] = vertices[3*(indicesInfo[3*(i+2)]-1)+j];
                  meshObject->_vertices[indicesInfo[3*(i+2)]-1]._normal[j] = vertexNorm[3*(indicesInfo[3*(i+2)+2]-1)+j];
                  meshObject->_vertices[indicesInfo[3*(i+2)]-1]._tangent[j] = vertexTan[3*(indicesInfo[3*(i+2)]-1)+j];
                  meshObject->_vertices[indicesInfo[3*(i+2)]-1]._bitangent[j] = vertexBiTan[3*(indicesInfo[3*(i+2)]-1)+j];
                }
              meshObject->_vertices[indicesInfo[3*(i+2)]-1]._uvCoord[0] = texCoord[2*(indicesInfo[3*(i+2)+1]-1)];
              meshObject->_vertices[indicesInfo[3*(i+2)]-1]._uvCoord[1] = texCoord[2*(indicesInfo[3*(i+2)+1]-1)+1];
	      
            }
	      
          numFacesRead += numIndicesInFace - 2;
          readingFaces = true;
        }
    }

  //Resize faces pointer as it could be mush larger than the actual number of faces
  meshObject->_numFaces = numFacesRead;
  meshObject->_faces = realloc(meshObject->_faces, sizeof(face) * numFacesRead);

  //Free all of the temporary arrays to store vertex data
  free(vertices);
  free(texCoord);
  free(vertexNorm);
  free(vertexTan);
  free(vertexBiTan);
  
  printf("Mesh %s loaded successfully\n",filename);

  //CalcMeshObjectVertexVectors(meshObject, true);
  
  return meshObject;
}

void NormalizeMeshObject(MeshObject *mo)
{
  assert(mo->_numVertices > 0);
  
  int i;

  float minX,minY,minZ,maxX,maxY,maxZ;

  minX = mo->_vertices[0]._xyzCoord[0]; maxX = mo->_vertices[0]._xyzCoord[0];
  minY = mo->_vertices[0]._xyzCoord[1]; maxY = mo->_vertices[0]._xyzCoord[1];
  minZ = mo->_vertices[0]._xyzCoord[2]; maxZ = mo->_vertices[0]._xyzCoord[2];

  //Find the min/max for x,y,z
  for (i=0; i<mo->_numVertices; i++)
    {
      if (mo->_vertices[i]._xyzCoord[0] < minX)
        {
          minX = mo->_vertices[i]._xyzCoord[0];
        }
      if (mo->_vertices[i]._xyzCoord[0] > maxX)
        {
          maxX = mo->_vertices[i]._xyzCoord[0];
        }
      if (mo->_vertices[i]._xyzCoord[1] < minY)
        {
          minY = mo->_vertices[i]._xyzCoord[1];
        }
      if (mo->_vertices[i]._xyzCoord[1] > maxY)
        {
          maxY = mo->_vertices[i]._xyzCoord[1];
        }
      if (mo->_vertices[i]._xyzCoord[2] < minZ)
        {
          minZ = mo->_vertices[i]._xyzCoord[2];
        }
      if (mo->_vertices[i]._xyzCoord[2] > maxZ)
        {
          maxZ = mo->_vertices[i]._xyzCoord[2];
        }
    }

  float scaleFactor = 1 / (max(maxX - minZ, max(maxY - minY, maxZ - minZ)));
  //Uniform scaling such that each point lies within 0 - 1
  for (i=0; i<mo->_numVertices; i++)
    {
      mo->_vertices[i]._xyzCoord[0] = scaleFactor*(mo->_vertices[i]._xyzCoord[0] - minX);
      mo->_vertices[i]._xyzCoord[1] = scaleFactor*(mo->_vertices[i]._xyzCoord[1] - minY);
      mo->_vertices[i]._xyzCoord[2] = scaleFactor*(mo->_vertices[i]._xyzCoord[2] - minZ);
    }
}

void CalcMeshObjectVertexVectors(MeshObject *mo, bool flipSign)
{
  int i,j;
  face *currFace;
  vec3_t vcross = vec3_create(NULL);
  vec3_t v1 = vec3_create(NULL);
  vec3_t v2 = vec3_create(NULL);

  //Zero out all the normals tangents and bitangents
  for (i=0; i<(mo->_numVertices); i++)
    {
      mo->_vertices[i]._normal[0] = 0;
      mo->_vertices[i]._normal[1] = 0;
      mo->_vertices[i]._normal[2] = 0;
    }

  // Find the normal of each face , add it to each vertex adjacent
  // To the face.
  for (i=0; i< mo->_numFaces; i++)
    {
      currFace = &(mo->_faces[i]);

      // Finding two vectors that represent the plane
      vec3_subtract(mo->_vertices[currFace->_indices[2]]._xyzCoord, mo->_vertices[currFace->_indices[0]]._xyzCoord, v1);
      vec3_subtract(mo->_vertices[currFace->_indices[1]]._xyzCoord, mo->_vertices[currFace->_indices[0]]._xyzCoord, v2);
      
      /* find cross-product between these vectors */
      vec3_cross(v1, v2, vcross);

      /* normalize this vector */
      vec3_normalize(v1,NULL);
      vec3_normalize(v2,NULL);
      vec3_normalize(vcross,NULL);

      /* add this normal to each vertex that is adjacent to face */
      for (j = 0; j < NUMVERTPERFACE; j++) {
        vec3_add(mo->_vertices[currFace->_indices[j]]._normal, vcross, NULL);
        vec3_add(mo->_vertices[currFace->_indices[j]]._tangent, v1, NULL);
        vec3_add(mo->_vertices[currFace->_indices[j]]._bitangent, v2, NULL);
      }
    }

  /* normalize all the normals at the vertices */

  for (i = 0; i < mo->_numVertices; i++)
    {
      if (flipSign)
        {
          vec3_negate(mo->_vertices[i]._normal, NULL);
          vec3_negate(mo->_vertices[i]._tangent, NULL);
          vec3_negate(mo->_vertices[i]._bitangent, NULL);
        }
      vec3_normalize(mo->_vertices[i]._normal, NULL);
      vec3_normalize(mo->_vertices[i]._tangent, NULL);
      vec3_normalize(mo->_vertices[i]._bitangent, NULL);
    }
}

void CalcPLYTextureCoordXY(MeshObject *mo)
{
  int i;
  for (i=0; i<mo->_numVertices; i++)
    {
      mo->_vertices[i]._uvCoord[0] = mo->_vertices[i]._xyzCoord[0];
      mo->_vertices[i]._uvCoord[1] = mo->_vertices[i]._xyzCoord[1];
    }
}

void CalcPLYTextureCoordYZ(MeshObject *mo)
{
  int i;
  for (i=0; i<mo->_numVertices; i++)
    {
      mo->_vertices[i]._uvCoord[0] = mo->_vertices[i]._xyzCoord[1];
      mo->_vertices[i]._uvCoord[1] = mo->_vertices[i]._xyzCoord[2];
    }
}

void CalcPLYTextureCoordXZ(MeshObject *mo)
{
  int i;
  for (i=0; i<mo->_numVertices; i++)
    {
      mo->_vertices[i]._uvCoord[0] = mo->_vertices[i]._xyzCoord[0];
      mo->_vertices[i]._uvCoord[1] = mo->_vertices[i]._xyzCoord[2];
    }
}

void SetMeshObjectMaterial(char * materialName, MeshObject *mo)
{
  Material *mat = getMaterial(materialName);

  if (mat == NULL)
    {
      fprintf(stderr, "Error: Material not set : %s could not be found\n",materialName);
      return;
    }

  mo->_material = mat;
}

/*
 *  Sends MeshObject vertices, normals, and faces
 *  to the GPU. Returns integer values specifying
 *  where the buffers for each one is located
 *  this allows for multiple objects to be created
 */
void initializeMeshObjectVertex(MeshObject *mo)
{
  glGenBuffers(1, &(mo->_vertexDataLoc));
  glGenBuffers(1, &(mo->_indexDataLoc));
  
  //Send vertex data to the GPU these include xyz, uv, normal, tangent, bitangent
  glBindBuffer(GL_ARRAY_BUFFER,mo->_vertexDataLoc);
  glBufferData(GL_ARRAY_BUFFER,sizeof(vertex)*mo->_numVertices,mo->_vertices,GL_STATIC_DRAW);
  
  //Setup the indicies for the faces to the GPU
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mo->_indexDataLoc);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(face)*mo->_numFaces,mo->_faces,GL_STATIC_DRAW);
  
  checkGLErrors("Loading ply mesh");
}

/*
 *  There are several requirements for drawing these objects in the shader
 *  If it is not drawing properly this is due to the shader not being complex
 *  enough.
 *  Vertex Properties: xyzcoor, uvcoor, normal, tangent, bitangent
 *  Textures: colormap, normalmap, displacementmap
 */
void DrawMeshObject(MeshObject *mo, shader *s)
{
  glUseProgram(s->_program);
  
  // Enable the verticies attributes array
  glBindBuffer(GL_ARRAY_BUFFER, mo->_vertexDataLoc);
  
  // Not sure if this need to be enabled each time
  // This will be moved if not nessisary
  // The shader program is smart enough to not
  // allocate an attribute index if value is not used
  //printf("%d %d %d %d %d\n",s->_vertexAttr[0],s->_vertexAttr[1],s->_vertexAttr[2],s->_vertexAttr[3],s->_vertexAttr[4]);
  glEnableVertexAttribArray(s->_vertexAttr[0]);
  glEnableVertexAttribArray(s->_vertexAttr[1]);
  glEnableVertexAttribArray(s->_vertexAttr[2]);
  glEnableVertexAttribArray(s->_vertexAttr[3]);
  glEnableVertexAttribArray(s->_vertexAttr[4]);
  
  // Locations are based upon stucture of vertex
  glVertexAttribPointer(s->_vertexAttr[0], 3, GL_FLOAT, GL_FALSE, sizeof(vertex), BUFFER_OFFSET(0));
  glVertexAttribPointer(s->_vertexAttr[1], 2, GL_FLOAT, GL_FALSE, sizeof(vertex), BUFFER_OFFSET(3 * sizeof(float)));
  glVertexAttribPointer(s->_vertexAttr[2], 3, GL_FLOAT, GL_TRUE, sizeof(vertex), BUFFER_OFFSET(5 * sizeof(float)));
  glVertexAttribPointer(s->_vertexAttr[3], 3, GL_FLOAT, GL_TRUE, sizeof(vertex), BUFFER_OFFSET(8 * sizeof(float)));
  glVertexAttribPointer(s->_vertexAttr[4], 3, GL_FLOAT, GL_TRUE, sizeof(vertex), BUFFER_OFFSET(11 * sizeof(float)));

  setOGLMaterial(mo->_material, s, GL_FRONT);
  
  // Bind the index array
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mo->_indexDataLoc);
  
  glDrawElements(GL_TRIANGLES, mo->_numFaces*NUMVERTPERFACE, GL_UNSIGNED_INT, NULL);

  checkGLErrors("End of plyobject draw function");
}

/*
 *  Prints the data contained within the PlyObject
 *  Use for relativly small ply files otherwise
 *  this will be HUGE. Make sure to call after normals
 *  have already been formed or it will SEGFAULT
 */
void PrintMeshObject(MeshObject *mo)
{
  int i;
  char buffer[128];
  printf("The Vertices: %d\n",mo->_numVertices);
  for (i=0; i<mo->_numVertices; i++)
    {
      printf("-Vertex %d:\n",i);
      vec3_str(mo->_vertices[i]._xyzCoord, buffer);
      printf("   XYZ Coordinates: %s\n",buffer);
      vec2_str(mo->_vertices[i]._uvCoord, buffer);
      printf("   UV Coordinates: %s\n",buffer);
      vec3_str(mo->_vertices[i]._normal, buffer);
      printf("   Normal Vector: %s\n",buffer);
      vec3_str(mo->_vertices[i]._tangent, buffer);
      printf("   Tangent Vector: %s\n",buffer);
      vec3_str(mo->_vertices[i]._bitangent, buffer);
      printf("   Bitangent Vector: %s\n",buffer);
    }
    
  printf("Faces: %d\n",mo->_numFaces);
  for (i=0; i<mo->_numFaces; i++)
    {
      
      printf("-Face %d: [%d %d %d]\n",i,mo->_faces[i]._indices[0],mo->_faces[i]._indices[1],mo->_faces[i]._indices[2]);
    }

  printf("Vertex and Face location:\n");
  printf("Vertex buffer location: %d\n", mo->_vertexDataLoc);
  printf("Face index buffer location: %d\n", mo->_indexDataLoc);

  printf("Material Properties:\n");
  if (mo->_material == NULL)
    {
      printf("No material properties defined for mesh\n");
    }
  else
    {
      printMaterial(mo->_material);
    }
}
 
