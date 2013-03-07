#ifndef TEXTURE_H
#define TEXTURE_H

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

typedef struct
{
  int _type;
  int _width;
  int _height;
  int _maxPixValue;
  char *_data;
}Texture;

Texture *readPPMFile(char *filename);
void writePPMFile(Texture *outputTexture, char *filename);
GLuint *initializeOGLTexture(Texture *t);

#endif
