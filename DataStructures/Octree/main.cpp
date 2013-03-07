#include "Octree.h"
#include "Point3d.h"

#include<time.h>
#include <iostream>
#include <stdlib.h>

int main(void)
{
  Octree octree;
  
  for (int i=0; i< 10; ++i)
    {
      Point3d p((1.0 * rand())/RAND_MAX,(1.0 * rand())/RAND_MAX,(1.0 * rand())/RAND_MAX);
      octree.Insert(p);
    }
  
  octree.PrintOctree();
  return 0;
}
