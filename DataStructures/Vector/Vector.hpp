#ifndef VECTOR_HPP
#define VECTOR_HPP

template <class T>
class Vector
{
public:
  int GetSize();
  void Insert(int index, T data);
  T At(int index);
  
protected:
  T *_dataVector;
  int _sizeVector;
}

#include "Vector.ipp"

#endif
