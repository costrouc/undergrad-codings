#ifndef SLNODE_HPP
#define SLNODE_HPP

#include <stdio.h>

template<class T>
class SLNode
{
 protected:
  SLNode<T> *_nextNode;
  T _dataValue;

 public:
  //Constructors
  SLNode();
  SLNode(SLNode<T> *nextNode);
  SLNode(T data);
  SLNode(SLNode<T> *nextNode, T data);

  // Getters
  T GetData();
  SLNode<T> *GetNext();
  
  // Setters
  void SetData(T data);
  void SetNext(SLNode<T> *nextNode);
};

#include "SLNode.ipp"

#endif
