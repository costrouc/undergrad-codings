#ifndef CUBE_H
#define CUBE_H

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

GLfloat cubeVerticies[] = {
  0.0,0.0,0.0,
  0.0,1.0,0.0,
  1.0,1.0,0.0,
  1.0,0.0,0.0,
  0.0,0.0,1.0,
  0.0,1.0,1.0,
  1.0,1.0,1.0,
  1.0,0.0,1.0
};

GLfloat cubeNormals[8][3] = {
  0,0,1,
  0,0,-1,
  1,0,0,
  -1,0,0,
  0,1,0,
  0,-1,0,
  0,-1,0,
  0,1,0
};

GLubyte cubeFaces[] = {
  4,7,6,5, 
  0,3,2,1, 
  3,2,6,7, 
  1,0,4,5, 
  5,6,2,1, 
  3,7,4,0, 
};

#endif
