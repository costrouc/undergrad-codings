#ifndef BSTNODE_HPP
#define BSTNODE_HPP

#include "DLNode.hpp"

template <class T>
class BSTNode : public DLNode<T>
{
 public:
  //Constructors
  BSTNode();
  BSTNode(BSTNode<T> *leftNode, BSTNode<T> *rightNode);
  BSTNode(T d);
  BSTNode(T d, BSTNode<T> *leftNode, BSTNode<T> *rightNode);

  //Getter Functions
  BSTNode<T> *GetLeft();
  BSTNode<T> *GetRight();

  //Setter Functions
  void SetLeft(BSTNode<T> *leftNode);
  void SetRight(BSTNode<T> *rightNode);
};

#include "BSTNode.ipp"

#endif
