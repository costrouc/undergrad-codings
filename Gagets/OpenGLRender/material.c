#include "material.h"
#include "texture.h"
#include "shader.h"
#include "utils.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

//Where the material library is stored, can only be seen in this file 
static Material *MaterialLibrary[10];
static int NumOfMaterials = 0;
static int MaxNumOfMaterials = 10;

Material *readMaterialFile(char *filename){
  Material *matObject = malloc(sizeof(Material));
  Texture *tempTexture;
  
  FILE *inputfile = fopen(filename, "r");
  if (inputfile == NULL)
    {
      fprintf(stderr,"Error: Failed to open material file: %s\n",filename);
    }

  printf("Attempting to Load material from .mtl file: %s\n",filename);
  
  char *lineBuffer = malloc(sizeof(char) * MAXLINELENGTH);
  char *tokenBuffer = malloc(sizeof(char) * MAXTOKENLENGTH);

  int i;
  
  while(fgets(lineBuffer, MAXLINELENGTH, inputfile) != NULL){
    tokenBuffer = strtok(lineBuffer, " \t\n");
    
    if(strcmp(tokenBuffer, "newmtl") == 0)
      {
        tokenBuffer = strtok(NULL, " \t\n");
        matObject->_name = malloc(sizeof(char) * strlen(tokenBuffer));
        strcpy(matObject->_name, tokenBuffer);
      }
    else if(strcmp(tokenBuffer, "Ka") == 0)
      {
        for(i=0; i<3; ++i)
          {
            tokenBuffer = strtok(NULL, " \t\n");
            matObject->_ambient[i] = atof(tokenBuffer);
          }
        matObject->_ambient[3] = 1.0;
      }
    else if(strcmp(tokenBuffer, "Kd") == 0)
      {
        for(i=0; i<3; ++i)
          {
            tokenBuffer = strtok(NULL, " \t\n");
            matObject->_diffuse[i] = atof(tokenBuffer);
          }
        matObject->_diffuse[3] = 1.0;
      }
    if(strcmp(tokenBuffer, "Ks") == 0)
      {
        for(i=0; i<3; ++i)
          {
            tokenBuffer = strtok(NULL, " \t\n");
            matObject->_specular[i] = atof(tokenBuffer);
          }
        matObject->_specular[3] = 1.0;
      }
    else if(strcmp(tokenBuffer, "Ns") == 0)
      {
        tokenBuffer = strtok(NULL, " \t\n");
        matObject->_shininess[0] = atof(tokenBuffer);
      }
    else if(strcmp(tokenBuffer, "map_Kd") == 0)
      {
        tokenBuffer = strtok(NULL, " \t\n");
        tempTexture = readPPMFile(tokenBuffer);
        if (tempTexture != NULL)
          {
            matObject->_texDataLoc[0] = initializeOGLTexture(tempTexture);
          }
      }
    else if(strcmp(tokenBuffer, "map_normal") == 0)
      {
        tokenBuffer = strtok(NULL, " \t\n");
        tempTexture = readPPMFile(tokenBuffer);
        if (tempTexture != NULL)
          {
            matObject->_texDataLoc[1] = initializeOGLTexture(tempTexture);
          }
      }
    else if(strcmp(tokenBuffer, "map_bump") == 0)
      {
        tokenBuffer = strtok(NULL, " \t\n");
        tempTexture = readPPMFile(tokenBuffer);
        if (tempTexture != NULL)
          {
            matObject->_texDataLoc[2] = initializeOGLTexture(tempTexture);
          }
        free(tempTexture);
      }
    else if(strcmp(tokenBuffer, "map_env_diff") == 0)
      {
        tokenBuffer = strtok(NULL, " \t\n");
        tempTexture = readPPMFile(tokenBuffer);
        if (tempTexture != NULL)
          {
            matObject->_envDataLoc[0] = initializeOGLTexture(tempTexture);
          }
        free(tempTexture);
      }
    else if(strcmp(tokenBuffer, "map_env_spec") == 0)
      {
        tokenBuffer = strtok(NULL, " \t\n");
        tempTexture = readPPMFile(tokenBuffer);
        if (tempTexture != NULL)
          {
            matObject->_envDataLoc[0] = initializeOGLTexture(tempTexture);
          }
        free(tempTexture);
      }
  }

  printf("Material %s loaded successfully\n",filename);
  
  return matObject;
}

void addMaterialToLibrary(Material *m)
{
  //Check if max number of materials have been added
  if (NumOfMaterials == MaxNumOfMaterials)
    {
      fprintf(stderr,"Error: Material was not added - Maximum number of materials is  %d\n", MaxNumOfMaterials);
      return;
    }

  //Check if material has previously been added
  int i;
  for (i=0; i<NumOfMaterials; i++)
    {
      if (strcmp(MaterialLibrary[i]->_name, m->_name) == 0)
        {
          fprintf(stderr,"Error: Material was not added - Material name %s already exits\n", m->_name);
          return;
        }
    }

  printf("Material %s added to material library\n", m->_name);

  MaterialLibrary[NumOfMaterials] = m;
  NumOfMaterials++;
}

Material *getMaterial(char *materialName)
{
  int i;
  for (i=0; i<NumOfMaterials; i++)
    {
      if (strcmp(MaterialLibrary[i]->_name, materialName) == 0)
        {
          return MaterialLibrary[i];
        }
    }
  
  fprintf(stderr,"Error: Material %s not found\n",materialName);
  return NULL;
}

void setOGLMaterial(Material *m, shader *s, GLenum face)
{
  // Enable the material properties
  glMaterialfv(face, GL_AMBIENT, m->_ambient);
  glMaterialfv(face, GL_DIFFUSE, m->_diffuse);
  glMaterialfv(face, GL_SPECULAR, m->_specular);
  glMaterialfv(face, GL_SHININESS, m->_shininess);

  // Define which texture loaction associates with each texture attribute
  glUniform1i(s->_texAttr[0], 0);
  glUniform1i(s->_texAttr[1], 1);
  glUniform1i(s->_texAttr[2], 2);

  // Enable the textures for the material
  glEnable(GL_TEXTURE_2D);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D,m->_texDataLoc[0]);
  glActiveTexture(GL_TEXTURE1);  
  glBindTexture(GL_TEXTURE_2D,m->_texDataLoc[1]);
  glActiveTexture(GL_TEXTURE2); 
  glBindTexture(GL_TEXTURE_2D,m->_texDataLoc[2]);

  // Define the enviroment textures
  glUniform1i(s->_envAttr[0], 3);
  glUniform1i(s->_envAttr[1], 4);

  glActiveTexture(GL_TEXTURE3);
  glBindTexture(GL_TEXTURE, m->_envDataLoc[0]);

  glActiveTexture(GL_TEXTURE4);
  glBindTexture(GL_TEXTURE, m->_envDataLoc[1]);
  
  checkGLErrors("End of setting opengl material properties for rendering");
}

void printMaterial(Material *m)
{
  printf("Material Name: %s\n",m->_name);
  printf("Ambient: %f %f %f %f\n",m->_ambient[0], m->_ambient[1], m->_ambient[2], m->_ambient[3]);
  printf("Diffuse: %f %f %f %f\n",m->_diffuse[0], m->_diffuse[1], m->_diffuse[2], m->_diffuse[3]);
  printf("Specular: %f %f %f %f\n",m->_specular[0], m->_specular[1], m->_specular[2], m->_specular[3]);
  printf("Shininess: %f\n", m->_shininess[0]);
  printf("Texture Data locations:\nDiffuse: %d\nNormal: %d\nBump %d\n", m->_texDataLoc[0],  m->_texDataLoc[1], m->_texDataLoc[2]);
  printf("EnvTexture Data locations:\nDiffuse: %d\nSpecular: %d\n", m->_envDataLoc[0], m->_envDataLoc[1]);
}
