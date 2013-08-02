#ifndef LIGHTS_H
#define LIGHTS_H

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

typedef struct
{
  char *_name;
  GLfloat _position[4];
  GLfloat _ambient[4];
  GLfloat _diffuse[4];
  GLfloat _specular[4];
}Light;

void setOGLLight(GLenum light,Light *l);

#endif
