#ifndef MATERIAL_H
#define MATERIAL_H

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

typedef struct
{
  char *_name;
  GLfloat _ambient[4];
  GLfloat _diffuse[4];
  GLfloat _specular[4];
  GLfloat _shininess[1];
}Material;

void setOGLMaterial(GLenum face, Material *m);

#endif
