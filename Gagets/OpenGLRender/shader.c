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

shader InitShader(char *vertFilename,char *fragFilename)
{
  char *vertShaderSource = readShaderFile(vertFilename);
  char *fragShaderSource = readShaderFile(fragFilename);
  
  GLuint vertShader;
  GLuint fragShader;
  shader s;
  
  s._program = glCreateProgram();
  
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
  
  glAttachShader(s._program, vertShader);
  glAttachShader(s._program, fragShader);

  glLinkProgram(s._program);
  glUseProgram(s._program);

  printf("Enabling: program\n");
  printProgramInfoLog(s._program);
  
  free(vertShaderSource);
  free(fragShaderSource);

  checkGLErrors("End of readShader()");
  return s;
}

void InitShaderAttrHandles(shader *s)
{
  glUseProgram(s->_program);
  
  s->_vertexAttr[0] = glGetAttribLocation(s->_program,"xyzcoord");
  s->_vertexAttr[1] = glGetAttribLocation(s->_program,"uvcoord");
  s->_vertexAttr[2] = glGetAttribLocation(s->_program,"normal");
  s->_vertexAttr[3] = glGetAttribLocation(s->_program,"tangent");
  s->_vertexAttr[4] = glGetAttribLocation(s->_program,"bitangent");

  s->_texAttr[0] = glGetUniformLocation(s->_program,"colormap");
  s->_texAttr[1] = glGetUniformLocation(s->_program,"normalmap");
  s->_texAttr[2] = glGetUniformLocation(s->_program,"bumpmap");

  s->_envAttr[0] = glGetUniformLocation(s->_program,"diffuse_irr_map");
  s->_envAttr[1] = glGetUniformLocation(s->_program,"specular_irr_map");
  
  checkGLErrors("End of shader hook loc of attr");
}

void printShader(shader *s)
{
  printf("Shader Information:\n");
  printf("Program: %d\n", s->_program);
  printf("Vertex Attribute locations:\n");
  printf("   Verticies: %d\n   Texture Coord: %d\n   Normals: %d\n   Tangents: %d\n   BiTangents %d\n",s->_vertexAttr[0], s->_vertexAttr[1], s->_vertexAttr[2], s->_vertexAttr[3], s->_vertexAttr[4]);
  printf("Texture Uniforms:\n");
  printf("   Diffuse Map: %d\n   Normal Map: %d\n   Bump Map: %d\n", s->_texAttr[0], s->_texAttr[1], s->_texAttr[2]);
  printf("Environment Texture Uniforms:\n");
  printf("   Diffuse Map: %d   Specular Map: %d\n", s->_envAttr[0], s->_envAttr[1]);
}

/*
 *  Prints an info log detailing the creation of a vertex
 *  of fragment shader. Used heavily for debuging purposes
 */
void printShaderInfoLog(GLuint s)
{
  GLint infoLogLength = 0;
  GLint charsWritten = 0;
  char *infoLog;

  glGetShaderiv(s, GL_INFO_LOG_LENGTH, &infoLogLength);

  if (infoLogLength > 0)
    {
      infoLog = malloc(sizeof(char) * infoLogLength);
      glGetShaderInfoLog(s, infoLogLength, &charsWritten, infoLog);
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

