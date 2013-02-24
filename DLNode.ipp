//Constructors
template <class T>
DLNode<T>::DLNode()
  : SLNode<T>(), _prevNode(NULL)                    
{
}

template <class T>
DLNode<T>::DLNode(DLNode<T> *prevNode, DLNode<T> *nextNode)
  : SLNode<T>(nextNode), _prevNode(prevNode)
{
}

template <class T>
DLNode<T>::DLNode(T d)
  : SLNode<T>(d),_prevNode(NULL)
{
}
  
template <class T> DLNode<T>::DLNode(T d, DLNode<T> *prevNode, DLNode<T> *nextNode)
  : SLNode<T>(d,nextNode), _prevNode(prevNode)
{
}

// Getters
template <class T>
DLNode<T> *DLNode<T>::GetPrev()
{
  return _prevNode;
}

// Setters
template <class T>
void DLNode<T>::SetPrev(DLNode<T> *prevNode)
{
  _prevNode = prevNode;
}
