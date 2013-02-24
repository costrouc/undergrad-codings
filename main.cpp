#include "BST.hpp"
#include "Queue.hpp"
#include "Stack.hpp"
#include "AVL.hpp"

#include <iostream>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  AVL<int> bst;

  for (int i=0; i<1000; i++)
    {
      int val = rand() % 1000;
      if (bst.insert(val) == true)
        {
          std::cout << "Inserted: " << val << std::endl;
        }
      else
        {
          std::cout << "Failed to Insert: " << val << std::endl;
        }
    }

  bst.preorder();
  bst.inorder();
  bst.postorder();
  bst.stats();
  
  for (int i=0; i<1000; i++)
    {
      int val = rand() % 100;
      if (bst.remove(val) == true)
        {
          std::cout << "Removed: " << val << std::endl;
        }
      else
        {
          std::cout << "Failed to Remove: " << val << std::endl;
        }
    }

  bst.preorder();
  bst.inorder();
  bst.postorder();
  
  Stack<int> s;
  s.Push(1);
  s.Push(2);
  s.Push(3);
  s.PrintStack();

  
  Queue<int> q;
  q.Enqueue(1);
  q.Enqueue(2);
  q.Enqueue(3);
  q.PrintQueue();

  
  return 0;
}
