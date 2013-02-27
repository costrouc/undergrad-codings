#ifndef TRIPLE_H
#define TRIPLE_H

#include <assert.h>

class Triple
{
 public:
  Triple();
  Triple(double x, double y, double z);
  ~Triple();
  
  double GetX();
  double GetY();
  double GetZ();
  void SetX(double x);
  void SetY(double y);
  void SetZ(double z);

  void Scale(double d);
  void Add(Triple t);
  bool Equals(Triple t); 
  bool EpsilonEquals(Triple t, double d);
  void ScaleAdd(double d,Triple t);
  void Absolute();
  
 protected:
  double _x;
  double _y;
  double _z;
};

#endif
