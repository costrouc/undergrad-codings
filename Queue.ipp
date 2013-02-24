template <class T>
Queue<T>::Queue()
  : _headNode(NULL), _tailNode(NULL), _numNodes(0)
{
}

template <class T>
void Queue<T>::Enqueue(T data)
{
  if (_headNode == NULL && _tailNode == NULL)
    {
      _headNode = _tailNode = new SLNode<T>(data);
    }
  else
    {
      SLNode<T> *tempNode = new SLNode<T>(data);
      _tailNode->SetNext(tempNode);
      _tailNode = tempNode;
    }
}

template <class T>
void Queue<T>::Dequeue()
{
  SLNode<T> *tempNode= _headNode;
  _headNode = _headNode->GetNext();
  delete(tempNode);
  _numNodes--;
}

template <class T>
T Queue<T>::Top()
{
  return _headNode;
}

template <class T>
int Queue<T>::GetSize()
{
  return _numNodes;
}

template <class T>
void Queue<T>::PrintQueue()
{
  SLNode<T> *traverseNode = _headNode;

  printf("Queue:\n");
  while(traverseNode != NULL)
    {
      std::cout << " " << traverseNode->GetData() << " ";
      traverseNode = traverseNode->GetNext();
    }
  std::cout << std::endl;
}
