/*
 * /brief Constructor for Stack().
 *
 * Initializes head node to NULL and size to 0.
 */
template <class T>
Stack<T>::Stack()
  : _headNode(NULL), _numNodes(0)
{
}

/*
 * /brief Destructor for Stack().
 *
 * Iteratively deletes each node by calling deconstructor for node. 
 */
template <class T>
Stack<T>::~Stack()
{
  SLNode<T> *nextNode    = _headNode;
  SLNode<T> *currentNode = _headNode;
  while (nextNode != NULL)
    {
      nextNode = currentNode->GetNext();
      delete(currentNode);
      currentNode = nextNode;
    }
}

/*
 * /brief Inserts elemenet of type T to begining of Stack.
 *
 * Allocates new node of type SLNode to begining.
 */
template <class T>
void Stack<T>::Push(T data)
{
  if (_headNode == NULL)
    {
      _headNode = new SLNode<T>(data);
    }
  else
    {
      SLNode<T> *tempNode = new SLNode<T>(_headNode,data);
      _headNode = tempNode;
    }
  
  _numNodes++;
}


template <class T>
T Stack<T>::Top()
{
  assert(!isEmpty());
  
  return _headNode->GetData();
}

template <class T>
void Stack<T>::Pop()
{
  assert(!isEmpty());
  
  SLNode<T> *tempNode = _headNode;
  _headNode = _headNode->GetNext();
  delete(tempNode);

  _numNodes--;
}

template <class T>
bool Stack<T>::isEmpty()
{
  return (_headNode == NULL);
}

template <class T>
int Stack<T>::GetSize()
{
  return _numNodes;
}

template <class T>
void Stack<T>::PrintStack()
{
  SLNode<T> *traverseNode = _headNode;

  printf("Stack:\n");
  while(traverseNode != NULL)
    {
      std::cout << " " << traverseNode->GetData() << " ";
      traverseNode = traverseNode->GetNext();
    }
  std::cout << std::endl;
}
