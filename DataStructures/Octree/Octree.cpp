#include "Octree.h"

Octree::Octree()
{  
  Point3d headNodeCenter(0.5,0.5,0.5);
  _headNode = OctreeNode(HEAD_NODE, 0, NULL, 0.5, headNodeCenter);
}

Octree::~Octree(){}

bool Octree::Insert(Point3d o)
{
  return _headNode.Insert(o);
}

std::vector<Point3d> Octree::GetObjects(Point3d p)
{
  return _headNode.GetObjects(p);
}

void Octree::PrintOctree()
{
  _headNode.PrintOctreeNode();
}
