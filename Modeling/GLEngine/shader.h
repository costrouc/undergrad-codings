#ifndef SHADER_H
#define SHADER_H

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif


char *read_shader_program(char *filename);
unsigned int InitShader(char *vertFile,char *fragFile);

void printShaderInfoLog(GLuint shader);
void printProgramInfoLog(GLuint shader);

#endif
