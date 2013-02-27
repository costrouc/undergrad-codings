// Constructor for Binary Search Tree (BST) 
template <class T>
BST<T>::BST()
{
  _rootNode = NULL;
  _numNodes = 0;
}

// insert - insert data value into BST. Base case BST is empty insert BSTNode
// to the root node. Otherwise use an itterative approach
template <class T>
bool BST<T>::insert(T data)
{
  if (isEmpty() == true)
    {
      _rootNode = new BSTNode<T>(data);
      _numNodes++;
      return true;
    }

  Stack<BSTNode<T>*> s;
  BSTNode<T> *currentNode = _rootNode;
  
  while (true)
    {
      s.Push(currentNode);
      int comparison = Comparator(data,currentNode->GetData());
      
      switch(comparison)
        {
        case LEFT:
          if (currentNode->GetLeft() == NULL)
            {
              currentNode->SetLeft(new BSTNode<T>(data));
              _numNodes++;
              Balance(&s);
              return true;
            }
          else
            {
              currentNode = currentNode->GetLeft();
            }
          break;
        case EQUAL:
          return false;
        case RIGHT:
          if (currentNode->GetRight() == NULL)
            {
              currentNode->SetRight(new BSTNode<T>(data));
              _numNodes++;
              Balance(&s);
              return true;
            }
          else
            {
              currentNode = currentNode->GetRight();
            }
          break;
        }
    }
}

template <class T>
bool BST<T>::remove(T data)
{
  if (isEmpty() == true)
    {
      return false;
    }
  
  void (BSTNode<T>::*setPrevNode)(BSTNode<T> *) = NULL;
  BSTNode<T> *previousNode = _rootNode;
  BSTNode<T> *currentNode = _rootNode;
  int comparison = Comparator(data,currentNode->GetData());;
  
  while (currentNode != NULL && comparison != EQUAL)
    {
      previousNode = currentNode;
      
      if (comparison == LEFT)
        {
          setPrevNode = &(BSTNode<T>::SetLeft);
          currentNode = currentNode->GetLeft();
        }
      else //==RIGHT
        {
          setPrevNode = &(BSTNode<T>::SetRight);
          currentNode = currentNode->GetRight();
        }
      if (currentNode == NULL)
        {
          break;
        }
      comparison = Comparator(data,currentNode->GetData());
    }
  
  if (currentNode == NULL)
    {
      //Element data not in BST
      return false;
    }
  else if (currentNode->GetRight() == NULL && currentNode->GetLeft() == NULL)
    {
      if (currentNode != _rootNode)
        {
          (previousNode->*setPrevNode)(NULL);
        }
      else
        {
          _rootNode = NULL;
        }
    }
  else if (currentNode->GetRight() == NULL)
    {
      if (currentNode != _rootNode)
        {
          (previousNode->*setPrevNode)(currentNode->GetLeft());
        }
      else
        {
          _rootNode = currentNode->GetLeft();
        }
    }
  else if (currentNode->GetLeft() == NULL)
    {
      if (currentNode != _rootNode)
        {
          (previousNode->*setPrevNode)(currentNode->GetRight());
        }
      else
        {
          _rootNode = currentNode->GetRight();
        }
    }
  else //Both left and right are not null
    {
      BSTNode<T> *currMaxLeftChild = currentNode->GetLeft();
      BSTNode<T> *prevMaxLeftChild = currentNode->GetLeft();

      //Get Max value of left tree
      while (currMaxLeftChild->GetRight() != NULL)
        {
          prevMaxLeftChild = currMaxLeftChild;
          currMaxLeftChild = currMaxLeftChild->GetRight();
        }

      currMaxLeftChild->SetRight(currentNode->GetRight());
      if (currMaxLeftChild != prevMaxLeftChild)
        {
          prevMaxLeftChild->SetRight(currMaxLeftChild->GetLeft());
          currMaxLeftChild->SetLeft(currentNode->GetLeft());
        }
      if (currentNode != _rootNode)
        {
          (previousNode->*setPrevNode)(currMaxLeftChild);
        }
      else
        {
          _rootNode = currMaxLeftChild;
        }
    }

  delete(currentNode);
  _numNodes--;
  return true;
}

template <class T>
void BST<T>::stats()
{
  int numNodes = 0;
  int numLeaves = 0;
  int numNodesOneChild = 0;
  int numEvenNodes = 0;
  int numOddNodes = 0;
  int sumHeightNodes = 0;
  int maxHeight = -1;
  int sumDepthNodes = 0;
  int maxDepth = -1;
  
  T sumNodeValues = 0;

  BSTNode<T> *currentNode;
  Stack<BSTNode<T>*> s;
  s.Push(_rootNode);
  
  while (s.isEmpty() == false)
    {
      currentNode = s.Top();

      //Check if a leaf or single parent
      if (currentNode->GetLeft() == NULL && currentNode->GetRight())
        {
          numLeaves++;
        }
      else if(currentNode->GetLeft() == NULL || currentNode->GetRight())
        {
          numNodesOneChild++;
        }

      //check if odd or even
      if (currentNode->GetData() % 2 == 0)
        {
          numEvenNodes++;
        }
      else
        {
          numOddNodes++;
        }

      sumNodeValues += currentNode->GetData();
      sumHeightNodes += GetHeight(currentNode);
      if (GetHeight(currentNode) > maxHeight)
        {
          maxHeight = GetHeight(currentNode);
        }

      if (GetDepth(currentNode) > maxDepth)
        {
          maxDepth = GetDepth(currentNode);
        }
      
      sumDepthNodes += GetDepth(currentNode);
      numNodes++;
      
      s.Pop();
      if (currentNode->GetRight() != NULL)
        {
          s.Push(currentNode->GetRight());
        }
      if (currentNode->GetLeft() != NULL)
        {
          s.Push(currentNode->GetLeft());
        }
    }

  printf("Tree Statistics:\n");
  printf("Number of Nodes in BST:        %d\n",numNodes);
  printf("Number of Leaves in BST:       %d\n",numLeaves);
  printf("Average height of the BST:     %.2f\n",(float) sumHeightNodes/(float) numNodes);
  printf("Max height of the BST:         %d\n",maxHeight);
  printf("Average depth of the BST:      %.2f\n",(float) sumDepthNodes/(float) numNodes);
  printf("Max depth of the BST:          %d\n",maxDepth);
  printf("Number of single parent nodes: %d\n",numNodesOneChild);
  printf("Number of even value nodes:    %d\n",numEvenNodes);
  printf("Number of odd value nodes:     %d\n",numOddNodes);
  printf("Sum nodes in BST:              ");
  std::cout << sumNodeValues << std::endl;
}

// Returns if an integer t is found within the BST
template <class T>
bool BST<T>::find(T data)
{
  BSTNode<T> *currentNode = _rootNode;
  int comparison;
  
  while (currentNode != NULL)
    {
      comparison = Comparator(data,currentNode->GetData());
      switch(comparison)
        {
        case LEFT:
          currentNode = currentNode->GetLeft();
          break;
        case EQUAL:
          return true;
        case RIGHT:
          currentNode = currentNode->GetRight();
          break;
        }          
    }
  return false;
}

// Returns NULL or the node in which an integer t is
// within the tree.
template <class T>
BSTNode<T> *BST<T>::findNode(T data)
{
  BSTNode<T> *currentNode = _rootNode;
  int comparison;
  
  while (currentNode != NULL)
    {
      comparison = Comparator(data,currentNode->GetData());
      switch(comparison)
        {
        case LEFT:
          currentNode = currentNode->GetLeft();
          break;
        case EQUAL:
          return currentNode;
        case RIGHT:
          currentNode = currentNode->GetRight();
          break;
        }          
    }
  
  return currentNode;
}

template <class T>
void BST<T>::printDotRep()
{
  Stack<BSTNode<T>*> s;
  BSTNode<T> *currentNode;
  s.Push(_rootNode);

  printf("digraph G {\n");
  while(s.isEmpty() == false)
    {
      currentNode = s.Top();
      s.Pop();
      if (currentNode->GetLeft() != NULL)
        {
          std::cout << currentNode->GetData() << " -> " << currentNode->GetLeft()->GetData() << std::endl;
          s.Push(currentNode->GetLeft());
        }
      if (currentNode->GetRight() != NULL)
        {
          std::cout << currentNode->GetData() << " -> " << currentNode->GetRight()->GetData() << std::endl;
          s.Push(currentNode->GetRight());
        }
    }
  printf("}\n");
}


// Returns the number of elements in the BST
template <class T>
int BST<T>::GetSize()
{
  return _numNodes;
}

// isEmpty - returns if tree is empty
template <class T>
bool BST<T>::isEmpty()
{
  return (_rootNode == NULL);
}

// Prints Tree in preorder NLR
template <class T>
void BST<T>::preorder()
{
  if (isEmpty() == true)
    {
      printf("Binary Search Tree (BST) is empty\n");
    }
  else
    {
      std::cout << "Elements in Preorder: ";
      Stack<BSTNode<T>*> s;
      BSTNode<T> *currentNode;
      s.Push(_rootNode);
      
      while (s.isEmpty() == false)
        {
          currentNode = s.Top();
          std::cout << " " << currentNode->GetData() << " ";
          s.Pop();
          if (currentNode->GetRight() != NULL)
            {
              s.Push(currentNode->GetRight());
            }
          if (currentNode->GetLeft() != NULL)
            {
              s.Push(currentNode->GetLeft());
            }
        }
      std::cout << std::endl;
    }
}

// Prints Tree in inorder LNR
template <class T>
void BST<T>::inorder()
{
  if (isEmpty() == true)
    {
      printf("Binary Search Tree (BST) is empty\n");
    }
  else
    {
      std::cout << "Elements in inorder: ";
      Stack<BSTNode<T>*> s;
      BSTNode<T> *currentNode = _rootNode;

      while (s.isEmpty() == false || currentNode != NULL)
        {
          if (currentNode != NULL)
            {
              s.Push(currentNode);
              currentNode = currentNode->GetLeft();
            }
          else
            {
              currentNode = s.Top();
              s.Pop();
              std::cout << " " << currentNode->GetData() << " ";
              currentNode = currentNode->GetRight();
            }
        }
      std::cout << std::endl;
    }
}

// Prints Tree in postorder LRN
template <class T>
void BST<T>::postorder()
{
  if (isEmpty() == true)
    {
      printf("Binary Search Tree (BST) is empty\n");
    }
  else
    {
      std::cout << "Elements in Postorder: ";
      Stack<BSTNode<T>*> s;
      Stack<BSTNode<T>*> output;
      BSTNode<T> *currentNode;
      s.Push(_rootNode);
      
      while (s.isEmpty() == false)
        {
          currentNode = s.Top();
          output.Push(currentNode);
          s.Pop();
          if (currentNode->GetLeft() != NULL)
            {
              s.Push(currentNode->GetLeft());
            }
          if (currentNode->GetRight() != NULL)
            {
              s.Push(currentNode->GetRight());
            }
        }
      while (output.isEmpty() == false)
        {
          std::cout << " " << output.Top()->GetData() << " ";
          output.Pop();
        }
      std::cout << std::endl;
    }
}

// Function used extensively to determine the position of an element within a tree
// Changing this funciton allows for the behavior of the BST to change easily.
template <class T> int BST<T>::Comparator(T data1, T data2)
{
  if (data1 > data2)
    {
      return RIGHT;
    }
  else if(data1 < data2)
    {
      return LEFT;
    }
  return EQUAL;
}

template <class T>
int BST<T>::GetHeight(BSTNode<T> *node)
{
  if (node == NULL)
    {
      return -1;
    }

  Stack<BSTNode<T>*> s;
  BSTNode<T>* currentNode;
  BSTNode<T>* previousNode = NULL;

  s.Push(node);
  
  int heightNode = -1;

  while(s.isEmpty() == false )
    {
      currentNode = s.Top();
      if (previousNode == NULL || previousNode->GetLeft() == currentNode || previousNode->GetRight() == currentNode)
        {
          if (currentNode->GetLeft() != NULL)
            {
              s.Push(currentNode->GetLeft());
            }
          else if (currentNode->GetRight() != NULL)
            {
              s.Push(currentNode->GetRight());
            }
        }
      else if(currentNode->GetLeft() == previousNode)
        {
          if (currentNode->GetRight() != NULL)
            {
              s.Push(currentNode->GetRight());
            }
        }
      else
        {
          s.Pop();
        }
      
      previousNode = currentNode;
      if (s.GetSize()-1 > heightNode)
        {
          heightNode = s.GetSize() - 1;
        }
    }
  return heightNode;
}

template <class T>
int BST<T>::GetDepth(BSTNode<T> *node)
{
  int currentDepth = -1;
  int comparison;
  BSTNode<T> *currentNode = _rootNode;

  while (currentNode != NULL && node != NULL)
    {
      currentDepth++;
      comparison = Comparator(node->GetData(),currentNode->GetData());
      switch(comparison)
        {
        case LEFT:
          currentNode = currentNode->GetLeft();
          break;
        case EQUAL:
          return currentDepth;
        case RIGHT:
          currentNode = currentNode->GetRight();
          break;
        }
    }
  return currentDepth;
}

template <class T>
void BST<T>::Balance(Stack<BSTNode<T>*> *s)
{
  printf("Regular BST's do not Balance the nodes!\n");
  printf("Parents: \n");
  while (s->isEmpty() == false)
    {
      std::cout << " " << s->Top()->GetData() << " ";
      s->Pop();
    }
  printf("\n");
}

