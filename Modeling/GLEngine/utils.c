#include "utils.h"

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <stdio.h>

void checkGLErrors(char *label)
{
  GLenum errorCode;
  const GLubyte *errorString;
  if ((errorCode = glGetError()) != GL_NO_ERROR)
  {
    errorString = gluErrorString(errorCode);
    printf("OpenGL Error: %s Label: %s\n",errorString,label);
  }
}
