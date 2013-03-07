#include "Ray.h"
Ray::Ray()
{
  _base = Point3d();
  _direction = Vector3d();
}

Ray::Ray(Point3d base, Vector3d direction)
{
  _base = base;
  _direction = direction;
}

Point3d Ray::GetBase()
{
  return _base;
}

Vector3d Ray::GetDirection()
{
  return _direction;
}
