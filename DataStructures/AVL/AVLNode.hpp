#ifndef AVLNODE_HPP
#define AVLNODE_HPP

#include "BSTNode.hpp"

template <class T>
class AVLNode : public BSTNode
{
public:
  //Constructors
  AVLNode();
  AVLNode(BSTNode<T> *leftNode, BSTNode<T> *rightNode);
  AVLNode(T d);
  AVLNode(T d, BSTNode<T> *leftNode, BSTNode<T> *rightNode);
  AVLNode(T d, BSTNode<T> *leftNode, BSTNode<T> *rightNode, int heightNode);
  
  //Getters
  int GetHeight();

  //Setters
  void SetHeight(int heightNode);

protected:
  int _heightNode;
};

#endif
