#ifndef STACK_HPP
#define STACK_HPP

#include "SLNode.hpp"

#include <iostream>
#include <assert.h>

template <class T>
class Stack
{
public:
  Stack();
  ~Stack();
  
  void Push(T data);
  void Pop();
  T Top();
  int GetSize();
  void PrintStack();
  bool isEmpty();
  
private:
  SLNode<T> *_headNode;
  int _numNodes;
};

#include "Stack.ipp"

#endif
