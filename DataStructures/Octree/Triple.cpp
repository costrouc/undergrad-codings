#include "Triple.h"

Triple::Triple()
{
  _x = 0;
  _y = 0;
  _z = 0;
}

Triple::Triple(double x, double y, double z)
{
  _x = x;
  _y = y;
  _z = z;
}

Triple::~Triple(){};

double Triple::GetX()
{
  return _x;
}

double Triple::GetY()
{
  return _y;
}

double Triple::GetZ()
{
  return _z;
}

void Triple::SetX(double x)
{
  _x = x;
}

void Triple::SetY(double y)
{
  _y = y;
}

void Triple::SetZ(double z)
{
  _z = z;
}

void Triple::Scale(double d)
{
  SetX(d*GetX());
  SetY(d*GetY());
  SetZ(d*GetZ());
}

void Triple::Add(Triple t)
{
  SetX(GetX() + t.GetX());
  SetY(GetY() + t.GetY());
  SetZ(GetZ() + t.GetZ());
}

bool Triple::Equals(Triple t)
{
  if (GetX() == t.GetX())
    {
      if (GetY() == t.GetY())
	{
	  if (GetZ() == t.GetZ())
	    {
	      return true;
	    }
	}
    }
  return false;
}

bool Triple::EpsilonEquals(Triple t, double d)
{
  assert(d >= 0);
  if (GetX() - t.GetX() < d && GetX() - t.GetX() > -d)
    {
      if (GetY() - t.GetY() < d && GetY() - t.GetY() > -d)
	{
	  if (GetZ() - t.GetZ() < d && GetZ() - t.GetZ() > -d)
	    {
	      return true;
	    }
	}
    }
  return false;
}

void Triple::ScaleAdd(double d, Triple t)
{
  Scale(d);
  Add(t);
}

void Triple::Absolute()
{
  if (GetX() < 0)
    {
      SetX(-GetX());
    }
  if (GetY() < 0)
    {
      SetY(-GetY());
    }
  if (GetZ() < 0)
    {
      SetZ(-GetZ());
    }
}
