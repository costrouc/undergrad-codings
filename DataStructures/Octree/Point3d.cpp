#include "Point3d.h"

Point3d::Point3d()
:Triple()
{
}

Point3d::Point3d(double x, double y, double z)
:Triple(x,y,z)
{
}

Point3d::~Point3d()
{
}


double Point3d::Distance(Point3d p)
{
  double dx = GetX()-p.GetX();
  double dy = GetY()-p.GetY();
  double dz = GetZ()-p.GetZ();
  
  return sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
}



