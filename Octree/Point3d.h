#ifndef POINT3D_H
#define POINT3D_H

#include "Triple.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>

class Point3d : public Triple
{
 public:
  Point3d();
  Point3d(double x,double y, double z);
  ~Point3d();
  
  double Distance(Point3d p);
};

#endif
