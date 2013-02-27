#include "Plane.h"

Plane::Plane()
{
  _point = Point3d();
  _normal = Vector3d();
}

Plane::Plane(Point3d point, Vector3d normal)
{
  _point = point;
  _normal = normal;
}

Plane::Plane(Plane &p)
{
  _point = p.GetPoint();
  _normal = p.GetNormal();
}

Point3d Plane::GetPoint()
{
  return _point;
}

Vector3d Plane::GetNormal()
{
  return _normal;
}

bool Plane::PointInPlane(Point3d p)
{
  Vector3d planeVec(p.GetX() - GetPoint().GetX(), p.GetY() - GetPoint().GetY(), p.GetZ() - GetPoint().GetZ());

  planeVec.Cross(_normal);

  Vector3d zeroVec(0,0,0);
  if (planeVec.EpsilonEquals(zeroVec,0.0000001))
    {
      return true;
    }
  return false;
}
