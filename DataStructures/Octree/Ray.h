#ifndef RAY_H
#define RAY_H

#include "Point3d.h"
#include "Vector3d.h"

/* Rays have not set length thus all rays have
 * a length of 1
 */

class Ray
{
  Ray();
  Ray(Point3d base, Vector3d direction);
  ~Ray();
  
  Point3d GetBase();
  Vector3d GetDirection();
  
 private:
  Point3d _base;
  Vector3d _direction;
};

#endif
