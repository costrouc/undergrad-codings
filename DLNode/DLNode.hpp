#ifndef DLNode_HPP
#define DLNode_HPP

#include "SLNode.hpp"

#include<stdio.h>

template <class T>
class DLNode : public SLNode<T>
{
 public:
  //Constructors
  DLNode();
  DLNode(DLNode<T> *prevNode, DLNode<T> *nextNode);
  DLNode(T d);
  DLNode(T d, DLNode<T> *prevNode, DLNode<T> *nextNode);

  //Getter Functions
  DLNode<T> *GetPrev();

  //Setter Functions
  void SetPrev(DLNode<T> *prevNode);
    
 protected:
  DLNode<T> *_prevNode;
};

#include "DLNode.ipp"

#endif
