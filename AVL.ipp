template <class T>
AVL<T>::AVL()
  :BST<T>()
{
}

template <class T>
void AVL<T>::Balance(Stack<BSTNode<T>*> *s)
{
  BSTNode<T> *currentNode;
  int leftHeight;
  int rightHeight;
  int balance;
  BSTNode<T> *retHeadNode = NULL;
  
  while (s->isEmpty() == false)
    {
      currentNode = s->Top();
      s->Pop();
      leftHeight = GetHeight(currentNode->GetLeft());
      rightHeight = GetHeight(currentNode->GetRight());
      balance = leftHeight - rightHeight;

      if (balance < -1)
        {
          if (GetHeight(currentNode->GetRight()->GetLeft()) - GetHeight(currentNode->GetRight()->GetRight()) > 0)
            {
              retHeadNode = rotateRightTwice(currentNode);
              break;
            }
          else
            {
              retHeadNode = rotateRightOnce(currentNode);
              break;
            } 
        }
      if (balance > 1)
        {
           if (GetHeight(currentNode->GetLeft()->GetLeft()) - GetHeight(currentNode->GetLeft()->GetRight()) < 0)
            {
              retHeadNode = rotateLeftTwice(currentNode);
              break;
            }
          else
            {
              retHeadNode = rotateLeftOnce(currentNode);
              break;
            } 
         }
     }
  
  if (retHeadNode != NULL)
    {
      if (s->isEmpty() == true)
        {
          this->_rootNode = retHeadNode;
        }
      else
        {
          BSTNode<T> *temp = s->Top();
          if (temp->GetRight() == currentNode)
            {
              temp->SetRight(retHeadNode);
            }
          else
            {
              temp->SetLeft(retHeadNode);
            }
        }
    }
  
   return;
}

template <class T>
BSTNode<T> *AVL<T>::rotateLeftOnce(BSTNode<T> *node1)
 {
   BSTNode<T> *node2 = node1->GetLeft();
   node1->SetLeft(node2->GetRight());
   node2->SetRight(node1);
   return node2;
}

template <class T>
BSTNode<T> *AVL<T>::rotateRightOnce(BSTNode<T> *node1)
{
  BSTNode<T> *node2 = node1->GetRight();
  node1->SetRight(node2->GetLeft());
  node2->SetLeft(node1);
  return node2;
}

template <class T>
BSTNode<T> *AVL<T>::rotateLeftTwice(BSTNode<T> *node1)
{
  node1->SetLeft(rotateRightOnce(node1->GetLeft()));
  return rotateLeftOnce(node1);
}

template <class T>
BSTNode<T> *AVL<T>::rotateRightTwice(BSTNode<T> *node1)
{
  node1->SetRight(rotateLeftOnce(node1->GetRight()));
  return rotateRightOnce(node1);
}
