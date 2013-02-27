// Constructors
template <class T>
SLNode<T>::SLNode()
  :_nextNode(NULL), _dataValue(NULL)
{
}

template <class T>
SLNode<T>::SLNode(SLNode<T> *nextNode)
  : _nextNode(nextNode), _dataValue(NULL)
{
}

template <class T>
SLNode<T>::SLNode(T data)
  : _nextNode(NULL), _dataValue(data)
{
}

template <class T>
SLNode<T>::SLNode(SLNode<T> *nextNode, T data)
  : _nextNode(nextNode), _dataValue(data)
{
}

//Getters
template <class T>
SLNode<T> *SLNode<T>::GetNext()
{
  return _nextNode;
}

template <class T>
T SLNode<T>::GetData()
{
  return _dataValue;
}

//Setters
template <class T>
void SLNode<T>::SetNext(SLNode<T> *nextNode)
{
  _nextNode = nextNode;
}

template <class T>
void SLNode<T>::SetData(T data)
{
  _dataValue = data;
}
