#ifndef BST_HPP
#define BST_HPP

#include "Stack.hpp"
#include "BSTNode.hpp"

#define RIGHT 1
#define EQUAL 0
#define LEFT -1

template <class T>
class BST
{
 public:
  BST();

  bool insert(T data);
  bool remove(T data);
  BSTNode<T> *findNode(T data);
  bool find(T data);
  bool isEmpty();
  void stats();
  
  void preorder();
  void inorder();
  void postorder();
  void printDotRep();
  
  int GetSize();
 protected:
  int Comparator(T data1, T data2);
  int GetHeight(BSTNode<T>* node);
  int GetDepth(BSTNode<T>* node);
  virtual void Balance(Stack<BSTNode<T>*> *s);
  
  BSTNode<T> *_rootNode;
  int _numNodes;
};

#include "BST.ipp"

#endif
