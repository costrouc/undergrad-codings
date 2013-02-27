#ifndef VECTOR3D_H
#define VECTOR3D_H

#include "Triple.h"

#include <math.h>

class Vector3d : public Triple
{
 public:
  Vector3d();
  Vector3d(double x, double y, double z);
  ~Vector3d();
  
  void Cross(Vector3d v);
  double Dot(Vector3d v);
  void Normalize();
  double Length();
};

#endif
