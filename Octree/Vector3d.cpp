#include "Vector3d.h"

Vector3d::Vector3d()
:Triple()
{
}

Vector3d::Vector3d(double x, double y, double z)
:Triple(x,y,z)
{
}

Vector3d::~Vector3d()
{
}

void Vector3d::Cross(Vector3d v)
{
  SetX(GetY() * v.GetZ() - GetX() * v.GetZ());
  SetY(GetZ() * v.GetX() - GetZ() * v.GetY());
  SetZ(GetX() * v.GetY() - GetY() * v.GetX());
}

double Vector3d::Dot(Vector3d v)
{
  double dx = GetX() * v.GetX();
  double dy = GetY() * v.GetY();
  double dz = GetZ() * v.GetZ();

  return dx+dy+dz;
}

void Vector3d::Normalize()
{
  Scale(1.0/Length());
}

double Vector3d::Length()
{
  double dx = pow(GetX(),2);
  double dy = pow(GetY(),2);
  double dz = pow(GetZ(),2);

  return sqrt(dx + dy + dz);
}

