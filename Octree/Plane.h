#ifndef PLANE_H
#define PLANE_H

#include "Point3d.h"
#include "Vector3d.h"

class Plane
{
 public:
  Plane();
  Plane(Point3d point, Vector3d normal);
  Plane(Plane &p);
  ~Plane();

  Point3d GetPoint();
  Vector3d GetNormal();

  bool PointInPlane(Point3d p);
  
 private:
  Point3d _point;
  Vector3d _normal;
};

#endif
