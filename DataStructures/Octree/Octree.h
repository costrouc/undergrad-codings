#ifndef OCTREE_H
#define OCTREE_H

#include "Point3d.h"
#include "OctreeNode.h"

#include <stdlib.h>
#include <iostream>
#include <vector>

/* Creates an octree that holds objects of type
 * Point3d. The box dimension are assumed to be (1,1,1).
 * A class Intersects is expected to have
 * A method IsIn() to determine if the object is
 * within the bounding box
 */

class Octree
{
 public:
  Octree();
  ~Octree();
  
  bool Insert(Point3d o);
  std::vector<Point3d> GetObjects(Point3d p);
  void PrintOctree();
  
 private:
  OctreeNode _headNode;
};

#endif
