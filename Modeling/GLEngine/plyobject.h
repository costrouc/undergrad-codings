#ifndef PLYOBJECT_H
#define PLYOBJECT_H

#include "material.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <stdbool.h>

#define MAXLINELENGTH 1024
#define MAXTOKENLENGTH 256

typedef struct
{
  GLfloat *_vertices;
  unsigned int _numVertices;
  GLfloat *_normals;
  unsigned int _numNormals;
  GLuint *_faces;
  unsigned int _numFaces;
  unsigned int _numVertPerFace;
  Material *_material;
}PLYObject;

PLYObject *readPLYFile(char *filename);
void CalcPLYObjectNormals(PLYObject *p, bool flipSign);
GLuint *initializePLYObject(PLYObject *p);
void NormalizePLYObject(PLYObject *p);
void DrawPLYObject(PLYObject *p);
void PrintPLYObject(PLYObject *p);

#endif
