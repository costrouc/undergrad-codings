//Constructors

template <class T>
AVLNode::AVLNode()
  :BSTNode(), _heightNode(NULL)
{
}

template <class T>
AVLNode::AVLNode(BSTNode<T> *leftNode, BSTNode<T> *rightNode)
  :BSTNode(leftNode,rightNode), _heightNode(NULL)
{
}

template <class T>
AVLNode::AVLNode(T data)
  :BSTNode(data), _heightNode(NULL)
{
}

template <class T>
AVLNode::AVLNode(T data, BSTNode<T> *leftNode, BSTNode<T> *rightNode)
  :BSTNode(data,leftNode,rightNode), _heightNode(NULL)
{
} 
  
template <class T>
AVLNode::AVLNode(T data, BSTNode<T> *leftNode, BSTNode<T> *rightNode, int heightNode)
  :BSTNode(data,leftNode,rightNode), _heightNode(heightNode)
{
}

//Getters

template <class T>
int AVLNode::GetHeight()
{
  return _heightNode;
}

//Setters

template <class T>
void AVLNode::SetHeight(int heightNode)
{
  _heightNode = heightNode;
}
