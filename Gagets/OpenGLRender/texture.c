#include "texture.h"
#include "utils.h"

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

static int MaxBufferSize = 1024;

/*
 *  Reads in a ppm file. Assumes that file is properly formated.
 *  May at some point add error checking.
 */
Texture *readPPMFile(char *filename)
{ 
  Texture *newTexture = malloc(sizeof(Texture));

  FILE *inputFile = fopen(filename,"r");

  if (inputFile == NULL)
    {
      fprintf(stderr,"Error: Failed to open input file %s\n", filename);
      free(newTexture);
      return NULL;
    }
  
  char *lineBuffer = malloc(sizeof(char) * MaxBufferSize);
  char *tokenBuffer = malloc(sizeof(char) * MaxBufferSize);

  bool readType = false;
  bool readDimensions = false;
  bool readMaxPixelValue = false; 

  printf("Attempting to read texture: %s\n",filename);
  
  while ( fgets(lineBuffer,MaxBufferSize,inputFile) != NULL)
    {
      tokenBuffer = strtok(lineBuffer," \t\n");
      
      if(readDimensions == true)
        {
          newTexture->_maxPixValue = atoi(tokenBuffer);
          readMaxPixelValue = true;

          int pixelsRead = fread(newTexture->_data, 3 * sizeof(char), newTexture->_width*newTexture->_height, inputFile);
          
          if (pixelsRead != newTexture->_width*newTexture->_height)
            {
              fprintf(stderr,"Error: was not able to read %d pixels only read %d\n", newTexture->_width*newTexture->_height, pixelsRead);
              free(newTexture);              
              return NULL;
            }
        }
      else if(readType == true)
        {
          if(strcmp(tokenBuffer,"#") != 0)
            {
              newTexture->_width = atoi(tokenBuffer);
              tokenBuffer = strtok(NULL," \t\n");
              newTexture->_height = atoi(tokenBuffer);

              newTexture->_data = malloc(sizeof(char) * 3 * newTexture->_width*newTexture->_height);
              readDimensions = true;
            }
        }
      else if (readType == false)
        {
          if (strcmp(tokenBuffer,"P6") == 0)
            {
              newTexture->_type = 6;
              readType = true;
            }
          else
            {
              fprintf(stderr,"Error: PPM in invalid format %s not supported yet\n",tokenBuffer);
              return NULL;
            }
        }
    }

  if (readMaxPixelValue != true)
    {
      fprintf(stderr,"Error: file %s is not in proper PPM format\n",filename);
      free(newTexture);
      return NULL;
    }

  printf("Texture %s loaded successfully\n",filename);
  
  return newTexture;
}

void writePPMFile(Texture *outputTexture, char *filename)
{
  FILE *outputFile = fopen(filename,"w");

  fprintf(outputFile,"P%d\n",outputTexture->_type);
  fprintf(outputFile,"# Created by Chris|Peter opengl software\n");
  fprintf(outputFile,"%d %d\n",outputTexture->_width,outputTexture->_height);
  fprintf(outputFile,"%d\n",outputTexture->_maxPixValue);
  fwrite(outputTexture->_data,3,outputTexture->_width*outputTexture->_height,outputFile);
}

GLuint initializeOGLTexture(Texture *t)
{
  GLuint bufTexture;
  glGenTextures(1,&bufTexture);

  //Copy texture to the buffer
  glBindTexture(GL_TEXTURE_2D,bufTexture);
  glTexImage2D(GL_TEXTURE_2D,0,GL_RGB,t->_width,t->_height,0,GL_RGB,GL_UNSIGNED_BYTE,t->_data);

  //Not sure what these do but if I want to adjust more parameters I can
  glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
  glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);

  checkGLErrors("End of Texture Initialization");
  
  //free(t) make sure to free the texture
  return bufTexture;
}
