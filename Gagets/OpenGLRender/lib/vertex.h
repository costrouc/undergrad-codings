#ifndef VERTEX_H
#define VERTEX_H

#include "glmathprototypes.h"

typedef struct
{
  float  _xyzCoord[3]; // Point coordinated
  float  _uvCoord[2];  // Texture Map Coordinates
  float _normal[3];	//vn
  float _tangent[3];	//vx
  float _bitangent[3];	//vy
}vertex;

#endif
