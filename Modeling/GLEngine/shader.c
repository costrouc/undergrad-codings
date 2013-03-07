#include "shader.h"
#include "utils.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/glext.h>
#else
#include <GL/glut.h>
#include <GL/glext.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>

char *readShaderFile(char *filename)
{
  FILE *inputFile;
  long fileLength;
  char *fileBuffer;
  
  inputFile = fopen(filename, "rb");

  if (inputFile == NULL)
    {
      fprintf(stderr,"Error: unable to open file: %s\n",filename);
      return NULL;
    }

  fseek(inputFile,0L, SEEK_END);
  fileLength = ftell(inputFile);
  rewind(inputFile);

  //Allocate memory to read in file
  fileBuffer = malloc(sizeof(char) * fileLength + 1);

  if (fileBuffer == NULL)
    {
      fclose(inputFile);
      fprintf(stderr,"Error: unable to allocate memory\n");
      return NULL;
    }

  if (1 != fread(fileBuffer, fileLength, sizeof(char), inputFile))
    {
      fclose(inputFile);
      fprintf(stderr,"Error: failed to read entire file %s\n",filename);
      return NULL;
    }

  fileBuffer[fileLength] = '\0';
  
  fclose(inputFile);
  return fileBuffer;
}

unsigned int InitShader(char *vertFilename,char *fragFilename)
{
  char *vertShaderSource = readShaderFile(vertFilename);
  char *fragShaderSource = readShaderFile(fragFilename);
  
  GLuint vertShader;
  GLuint fragShader;
  GLuint program;

  program = glCreateProgram();
  
  vertShader = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vertShader, 1, (const char**) &vertShaderSource, NULL);
  glCompileShader(vertShader);
  printf("Compiling: %s\n",vertFilename);
  printShaderInfoLog(vertShader);

  fragShader = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fragShader, 1, (const char**) &fragShaderSource, NULL);
  glCompileShader(fragShader);
  printf("Compiling: %s\n",fragFilename);
  printShaderInfoLog(fragShader);
  
  glAttachShader(program, vertShader);
  glAttachShader(program, fragShader);

  glLinkProgram(program);
  glUseProgram(program);

  printf("Enabling: program\n");
  printProgramInfoLog(program);
  //checkGLErrors("End of readShader()");
  free(vertShaderSource);
  free(fragShaderSource);
  return program;
}

/*
 *  Prints an info log detailing the creation of a vertex
 *  of fragment shader. Used heavily for debuging purposes
 */
void printShaderInfoLog(GLuint shader)
{
  GLint infoLogLength = 0;
  GLint charsWritten = 0;
  char *infoLog;

  glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLogLength);

  if (infoLogLength > 0)
    {
      infoLog = malloc(sizeof(char) * infoLogLength);
      glGetShaderInfoLog(shader, infoLogLength, &charsWritten, infoLog);
      printf("%s\n",infoLog);
      free(infoLog);
    }
}

/*
 *  Prints an info log detailing the creation of a GPU program
 *  (vertex + fragment). Used heavily for debuging purposes.
 */
void printProgramInfoLog(GLuint program)
{
  GLint infoLogLength = 0;
  GLint charsWritten = 0;
  char *infoLog;

checkGLErrors("End of readShader()");

  glGetProgramiv(program, GL_INFO_LOG_LENGTH, &infoLogLength);

  if (infoLogLength > 0)
    {
      infoLog = malloc(sizeof(char) * infoLogLength);
      glGetProgramInfoLog(program, infoLogLength, &charsWritten, infoLog);
      printf("%s\n",infoLog);
      free(infoLog);
    }
}

