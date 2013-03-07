#ifndef QUEUE_HPP
#define QUEUE_HPP

#include "SLNode.hpp"

template <class T>
class Queue
{
public:
  Queue();

  void Enqueue(T data);
  void Dequeue();
  T Top();
  int GetSize();
  void PrintQueue();
  
protected:
  SLNode<T> *_headNode;
  SLNode<T> *_tailNode;
  int _numNodes;
};

#include "Queue.ipp"

#endif
