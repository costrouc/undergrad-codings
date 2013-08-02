#ifndef SHADER_H
#define SHADER_H

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

typedef struct
{
  GLuint _program;
  GLuint _vertexAttr[5];
  GLuint _texAttr[3];
  GLuint _envAttr[2];
}shader;


char *read_shader_program(char *filename);
shader InitShader(char *vertFile,char *fragFile);
void InitShaderAttrHandles(shader *s);

void printShader(shader *s);
void printShaderInfoLog(GLuint shader);
void printProgramInfoLog(GLuint shader);





#endif
