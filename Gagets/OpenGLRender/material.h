#ifndef MATERIAL_H
#define MATERIAL_H

#include "shader.h"

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
  GLuint _texDataLoc[3];
  GLuint _envDataLoc[2];
}Material;

Material *readMaterialFile(char *filename);

void addMaterialToLibrary(Material *m);

Material *getMaterial(char *materialName);
void setOGLMaterial(Material *m, shader *s, GLenum face);
void printMaterial(Material *m);

#endif
