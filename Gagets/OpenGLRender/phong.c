// Main file to execute the ply rendering calls all of the other
// major functions

#include "utils.h"
#include "texture.h"
#include "mouseinput.h"
#include "keyinput.h"
#include "light.h"
#include "meshobject.h"
#include "shader.h"
#include "phong.h"
#include "computer.h"
#include "lib/glmath.h"

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
#include <math.h>

computer PCsetting = {
  ._windowHeight = 400,
  ._windowWidth = 400,
  ._windowOffsetX = 100,
  ._windowOffsetY = 100
};

Light keyLight = {
  ._name = "keyLight",
  ._position = {3.0, 3.0, 3.0, 1.0},
  ._ambient = {0.5, 0.5, 0.5, 1.0},
  ._diffuse = {1.0, 1.0, 1.0, 1.0},
  ._specular = {0.5, 0.5, 0.5, 1.0}
};

Light fillLight = {
  ._name = "fillLight",
  ._position = {-3.0, -3.0, 2.0, 1.0},
  ._ambient = {0.2, 0.2, 0.2, 1.0},
  ._diffuse = {1.0, 1.0, 1.0, 1.0},
  ._specular = {0.2, 0.2, 0.2, 1.0}
};

Light backLight = {
  ._name = "backLight",
  ._position = {-2.0, -2.0, 3.0, 1.0},
  ._ambient = {0.2, 0.2, 0.2, 1.0},
  ._diffuse = {1.0, 1.0, 1.0, 1.0},
  ._specular = {0.2, 0.2, 0.2, 1.0}
};

//The initial setup of the eye
float eye[] = {3.0,3.0,3.0};
float viewpt[] = {0.0,0.0,0.0};
float up[] = {0.0,1.0,0.0};

MeshObject *meshObject;
MeshObject *planeObject;

GLuint envTexture;

shader *shaderPrograms;

int moving = 1;
float objRot_x = 0.0;
float objRot_y = 0.0;
float objRot_z = 0.0;
float trans_x = 0.0;
float trans_y = 0.0;
float trans_z = 0.0;

void setupViewVolume()
{ 
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0,1.0,0.5,100.0);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(eye[0],eye[1],eye[2],viewpt[0],viewpt[1],viewpt[2],up[0],up[1],up[2]);
}

void DrawEnviromentTexture()
{
  GLUquadric *qptr = gluNewQuadric();

  glUseProgram(0);
  gluQuadricTexture(qptr, 1);
  gluQuadricOrientation(qptr, GLU_INSIDE);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, envTexture);
  glEnable(GL_TEXTURE_2D);
  gluSphere(qptr, 10.0, 64, 64);

  checkGLErrors("End of draw enviroment sphere");
}

void renderScene(void)
{  
  glClearColor(0.0, 0.0, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      
  int n = 1; //number of light rays to shoot
  float aperture = 0.03;

  vec3_t toEyeVec = vec3_create(NULL);
  vec3_subtract(viewpt, eye, toEyeVec);
  vec3_t right = vec3_create(NULL);
  vec3_cross(toEyeVec, up, right);
  vec3_normalize(right,NULL);
  vec3_t p_up = vec3_create(NULL);
  vec3_cross(toEyeVec, right, p_up);
  vec3_normalize(p_up, NULL);

  vec3_t component1 = vec3_create(NULL);
  vec3_t component2 = vec3_create(NULL);
  vec3_t bokeh = vec3_create(NULL);
  
  int i;
  for(i=0; i<n; i++)
    {
      vec3_scale(right, cosf(i * 2 * M_PI / n), component1);
      vec3_scale(p_up, sinf(i * 2 * M_PI / n), component1);
      vec3_add(component1, component2, bokeh);
      vec3_scale(bokeh, aperture, NULL);

      glLoadIdentity();
      gluLookAt(eye[0]+bokeh[0], eye[1]+bokeh[1], eye[2]+bokeh[2],viewpt[0],viewpt[1],viewpt[2],up[0], up[1], up[2]);
      
      /* Object Rotations/Translations*/
      glPushMatrix();

         glTranslatef(0.0 + trans_x, 0.0 + trans_y, 0.0 + trans_z);
	 glRotatef((360.0 / (30 * 1)) * objRot_x, 1, 0, 0);
	 glRotatef((360.0 / (30 * 1)) * objRot_y, 0, 1, 0);
	 glRotatef((360.0 / (30 * 1)) * objRot_z, 0, 0, 1);

	 glScalef(1.0, 1.0, 1.0);

	 glEnable(GL_LIGHTING);
	 glEnable(GL_LIGHT0);
	 glEnable(GL_LIGHT1);
	 glEnable(GL_LIGHT2);
	 
	 //Draw subject Object
	 DrawMeshObject(meshObject, &(shaderPrograms[2]));
	
	DrawMeshObject(planeObject, &(shaderPrograms[2]));

	 glDisable(GL_LIGHTING);
	 glDisable(GL_LIGHT0);
	 glDisable(GL_LIGHT1);
	 glDisable(GL_LIGHT2);
	 
	 //Draw Enviroment Texture
	 DrawEnviromentTexture();

	 glAccum(i ? GL_ACCUM : GL_LOAD, 1.0 / n);
      glPopMatrix();
    }
  glAccum(GL_RETURN, 1);
    
  glutSwapBuffers();

  checkGLErrors("Display Function");
}

void initOGL(int argc,char **argv)
{
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | GLUT_ACCUM);
  glutInitWindowPosition(PCsetting._windowOffsetX, PCsetting._windowOffsetY);
  glutInitWindowSize(PCsetting._windowWidth, PCsetting._windowHeight);
  glutCreateWindow("Start of Chris | Peter Project");
  glEnable(GL_DEPTH_TEST);

  // Set the starting view point
  setupViewVolume();

  //Enable and set all the lights
  setOGLLight(GL_LIGHT0,&keyLight);
  setOGLLight(GL_LIGHT1,&fillLight);
  setOGLLight(GL_LIGHT2,&backLight);

  Texture *envtex = readPPMFile("textures/space.ppm");
  envTexture = initializeOGLTexture(envtex);

  // Initialize materials library
  Material *tempMaterial = readMaterialFile("matfiles/plane.mtl");
  addMaterialToLibrary(tempMaterial);

  // Enable Ply Object
  planeObject = readPLYFile("plyfiles/airplane.ply");
  CalcPLYTextureCoordXZ(planeObject);
  //PrintMeshObject(meshObject);
  initializeMeshObjectVertex(planeObject);
  SetMeshObjectMaterial("look2", planeObject);  
  
  //Now lets test the obj reader and set up a plane
  meshObject = readPLYFile("plyfiles/airplane.ply");
  //NormalizeMeshObject(meshObject);
  //CalcMeshObjectVertexVectors(meshObject, true);
  initializeMeshObjectVertex(meshObject);
  //PrintMeshObject(meshObject);
  
  shaderPrograms = malloc(sizeof(shader) * 3);
  shaderPrograms[0]._program = 0;
  shaderPrograms[2] = InitShader("shader/texturevert.glsl","shader/texturefrag.glsl");
  InitShaderAttrHandles(&(shaderPrograms[2]));
  printShader(&(shaderPrograms[2]));
}
