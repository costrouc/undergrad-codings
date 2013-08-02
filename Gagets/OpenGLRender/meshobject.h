#ifndef MESHOBJECT_H
#define MESHOBJECT_H

#include "material.h"
#include "lib/vertex.h"
#include "lib/face.h"
#include "shader.h"
#include "utils.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <stdbool.h>

typedef struct
{
  int _numVertices;
  vertex *_vertices;
  int _numFaces;
  face *_faces;
  GLuint _vertexDataLoc;
  GLuint _indexDataLoc;
  Material *_material;
}MeshObject;

MeshObject *readPLYFile(char *filename);
MeshObject *readOBJFile(char *filename);

void CalcMeshObjectVertexVectors(MeshObject *p, bool flipSign);
void CalcPLYTextureCoordXY(MeshObject *p);
void CalcPLYTextureCoordYZ(MeshObject *p);
void CalcPLYTextureCoordXZ(MeshObject *p);

void SetMeshObjectMaterial(char *materialName, MeshObject *p);

void NormalizeMeshObject(MeshObject *p);

void initializeMeshObjectVertex(MeshObject *p);
void DrawMeshObject(MeshObject *p, shader *s);

void PrintMeshObject(MeshObject *p);

#endif
