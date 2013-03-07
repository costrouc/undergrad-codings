//
// Phong shading, with lighting done in eye coordinates, thater than
// world coordinates.
//

#include "utils.h"
#include "texture.h"
#include "mouseinput.h"
#include "keyinput.h"
#include "light.h"
#include "plyobject.h"
#include "shader.h"
#include "phong.h"


#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <OpenGL/glext.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>
#endif

#include <stdio.h>
#include <stdlib.h>

Light keyLight = {
    ._name = "keyLight",
    ._position = {0.0, 1.0, 1.75, 0.5},
    ._ambient = {0.6, 0.2, 1.0, 0.67},
    ._diffuse = {0.4, 0.0, 1.0, 0.5},
    ._specular = {0.9, 0.9, 0.9, 0.2}
  };

Light fillLight = {
  ._name = "fillLight",
  ._position = {1.0, 1.5, 0.35, 0.75},
  ._ambient = {0.2, 0.2, 1.0, 0.6},
  ._diffuse = {0.2, 0.4, 1.0, 0.5},
  ._specular = {0.9, 0.9, 0.9, 0.2}
};

Light backLight = {
  ._name = "backLight",
  ._position = {-2.0, 1.25, 0.0, 0.8},
  ._ambient = {0.8, 0.0, 1.0, 0.5},
  ._diffuse = {0.8, 0.0, 1.0, 0.5},
  ._specular = {0.9, 0.9, 0.9, 0.2}
};

Material silver = {
  ._name = "silver",
  ._ambient = {0.1, 0.1, 0.1, 1.0},
  ._diffuse = {1.0, 1.0, 1.0, 1.0},
  ._specular = {0.5, 0.5, 0.5, 1.0},
  ._shininess = {1.8}
};

PLYObject *plyObject;
// attrLocation 0 -> Verticies
// attrLocation 1 -> Normals
// attrLocation 2 -> Faces
GLuint *attrLocation;
// textureLocation 0 -> only texture
GLuint *textureLocation;
GLuint *textureNormalLocation;
GLuint *textureDispLocation;

GLuint shaderProgDiffNormDisp;
GLuint shaderProgPhong;
GLuint shaderProgNormal = 0;

int moving = 1;
float objRot_x = 0.0;
float objRot_y = 0.0;
float objRot_z = 0.0;
float trans_x = 0.0;
float trans_y = 0.0;
float trans_z = 0.0;

void setupViewVolume()
{
  //The initial setup of the eye
  GLfloat eye[] = {2.0,2.0,2.0};
  GLfloat viewpt[] = {0.0,0.0,0.0};
  GLfloat up[] = {0.0,1.0,0.0};
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0,1.0,0.5,20.0);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(eye[0],eye[1],eye[2],viewpt[0],viewpt[1],viewpt[2],up[0],up[1],up[2]);
}

void renderScene(void)
{
    
  glClearColor(0.5, 0.5, 0.5, 1.0);

  glUseProgram(shaderProgPhong);

  /*
   *  Camara movement should be put right here! I finally understand
   */

  glClearAccum(0.0,0.0,0.0,1.0);
  glClear(GL_ACCUM_BUFFER_BIT);

  int aliasingiter = 1;
  float scaleJitter = 0.0001;
  int i;
  for (i=0; i<aliasingiter; i++)
    {
      //Need to clear out color buffer but between accumulations
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      
      /* Object Rotations/Translations*/
      glPushMatrix();
      //Jitter the point for anti-aliasing must be here
      glTranslatef(scaleJitter*(float)rand()/(float)RAND_MAX, scaleJitter*(float)rand()/(float)RAND_MAX, scaleJitter*(float)rand()/(float)RAND_MAX);
      
      glTranslatef(0.0 + trans_x, 0.0 + trans_y, 0.0 + trans_z);
      glRotatef((360.0 / (30 * 1)) * objRot_x, 1, 0, 0);
      glRotatef((360.0 / (30 * 1)) * objRot_y, 0, 1, 0);
      glRotatef((360.0 / (30 * 1)) * objRot_z, 0, 0, 1);
      glScalef(1.0, 1.0, 1.0);
   
      if (plyObject->_numVertPerFace == 3)
        {
          glDrawElements(GL_TRIANGLES,plyObject->_numFaces*plyObject->_numVertPerFace,GL_UNSIGNED_INT,NULL);
        }
      else if (plyObject->_numVertPerFace == 4)
        {
          glDrawElements(GL_QUADS,plyObject->_numFaces*plyObject->_numVertPerFace,GL_UNSIGNED_INT,NULL);
        }

      glUseProgram(shaderProgDiffNormDisp);   
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D,textureLocation[0]);
      glActiveTexture(GL_TEXTURE1);  
      glBindTexture(GL_TEXTURE_2D,textureNormalLocation[0]);
      glActiveTexture(GL_TEXTURE2); 
      glBindTexture(GL_TEXTURE_2D,textureDispLocation[0]); 

      glEnable(GL_TEXTURE_2D);
      glBegin(GL_QUADS);
      glTexCoord2f(0.0,1.0);
      glVertex3f(-1.0,0.1775,-1.0);
      glTexCoord2f(1.0,1.0);
      glVertex3f(-1.0,0.1775,1.0);
      glTexCoord2f(1.0,0.0);
      glVertex3f(1.0,0.1775,1.0);
      glTexCoord2f(0.0,0.0);
      glVertex3f(1.0,0.1775,-1.0);
      glEnd();
      glDisable(GL_TEXTURE_2D);
   
      glPopMatrix();

      glAccum(GL_ACCUM, 1.f/aliasingiter);
    }

  glAccum(GL_RETURN, 1.f);
  glFlush();

  glutSwapBuffers();

  checkGLErrors("Display Function");
}

void initOGL(int argc,char **argv)
{
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | GLUT_ACCUM);
  glutInitWindowPosition(100,100);
  glutInitWindowSize(512,512);
  glutCreateWindow("Start of Chris | Peter Project");
  glEnable(GL_DEPTH_TEST);
  setupViewVolume();

  // Materials Predefined in light.h and material.h
  
  setOGLLight(GL_LIGHT0,&keyLight);
  setOGLLight(GL_LIGHT1,&fillLight);
  setOGLLight(GL_LIGHT2,&backLight);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glEnable(GL_LIGHT2); 

  setOGLMaterial(GL_FRONT,&silver);

  shaderProgPhong = InitShader("shader/phongvert.glsl", "shader/phongfrag.glsl");

  //Argv[1] is specified to be the filename
  plyObject = readPLYFile(argv[1]);
  NormalizePLYObject(plyObject);
  //PrintPLYObject(plyObject);
  attrLocation = initializePLYObject(plyObject);

  shaderProgDiffNormDisp = InitShader("shader/texturevert.glsl","shader/texturefrag.glsl");
  
  glUseProgram(shaderProgDiffNormDisp);
  //Get location of uniform variables in shader
  GLint textureLoc = glGetUniformLocation(shaderProgDiffNormDisp,"textureMap");
  GLint normalLoc = glGetUniformLocation(shaderProgDiffNormDisp,"normalMap");
  GLint dispLoc = glGetUniformLocation(shaderProgDiffNormDisp,"dispMap");
  glUniform1i(textureLoc, 0);
  glUniform1i(normalLoc, 1);
  glUniform1i(dispLoc, 2);

  //Create Textures
  Texture *texture1 = readPPMFile("textures/pebbles.ppm");
  //writePPMFile(texture,"textures/exampleoutput.ppm");
  textureLocation = initializeOGLTexture(texture1);

  Texture *texture2 = readPPMFile("textures/pebblesnormal.ppm");
  textureNormalLocation = initializeOGLTexture(texture2);

  Texture *texture3 = readPPMFile("textures/pebblesdisplmap.ppm");
  textureDispLocation = initializeOGLTexture(texture3);
}
