#ifndef AVL_HPP
#define AVL_HPP

#include "BST.hpp"
#include "BSTNode.hpp"
#include "Stack.hpp"

template <class T>
class AVL : public BST<T>
{
public:
  AVL(); 

protected:
  void  Balance(Stack<BSTNode<T>*> *s);
  BSTNode<T> *rotateLeftOnce(BSTNode<T> *node1);
  BSTNode<T> *rotateRightOnce(BSTNode<T> *node1);
  BSTNode<T> *rotateRightTwice(BSTNode<T> *node1);
  BSTNode<T> *rotateLeftTwice(BSTNode<T> *node1);
};

#include "AVL.ipp"

#endif
