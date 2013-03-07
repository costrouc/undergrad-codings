//Constructors

template <class T>
BSTNode<T>::BSTNode()
  : DLNode<T>()
{
}

template <class T>
BSTNode<T>::BSTNode(BSTNode<T> *leftNode, BSTNode<T> *rightNode)
  : DLNode<T>(leftNode,rightNode)
{
}

template <class T>
BSTNode<T>::BSTNode(T d)
  : DLNode<T>(d)
{
}

template <class T>
BSTNode<T>::BSTNode(T d, BSTNode<T> *leftNode, BSTNode<T> *rightNode)
  : DLNode<T>(d,leftNode,rightNode)
{
}

//Getter Functions

template <class T>
BSTNode<T> *BSTNode<T>::GetLeft()
{
  return (BSTNode<T> *) this->_prevNode;
}

template <class T>
BSTNode<T> *BSTNode<T>::GetRight()
{
  return (BSTNode<T> *) this->_nextNode;
}

//Setter Functions

template <class T>
void BSTNode<T>::SetLeft(BSTNode<T> *leftNode)
{
  this->_prevNode = leftNode;
}

template <class T>
void BSTNode<T>::SetRight(BSTNode<T> *rightNode)
{
  this->_nextNode = rightNode;
}
