#include "Sphere.h"

Sphere::Sphere()
{
  _center = Point3d();
  _radius = 0;
}

Sphere::Sphere(Point3d center, double radius)
{
  _center = center;
  _radius = radius;
}

Sphere::~Sphere()
{
}

Point3d Sphere::GetCenter()
{
  return _center;
}

double Sphere::GetRadius()
{
  return _radius;
}
