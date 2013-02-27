#ifndef SPHERE_H
#define SPHERE_H

#include "Point3d.h"

class Sphere
{
  Sphere();
  Sphere(Point3d center, double radius);
  ~Sphere();
  
  Point3d GetCenter();
  double GetRadius();

 private:
  Point3d _center;
  double _radius;
};

#endif
