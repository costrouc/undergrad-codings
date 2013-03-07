#ifndef OCTREENODE_H
#define OCTREENODE_H

#include <assert.h>
#include <vector>
#include <iostream>
#include <stdlib.h>

#include "Point3d.h"

#define NUM_OF_CHILDREN 8 //Do not change
#define MAX_DEPTH 7
#define MAX_OBJECTS_PER_NODE 8
#define MIN_CONTAINER_SIZE 0.000001

enum node_types_t {HEAD_NODE, UP_RIGHT_FRONT, UP_RIGHT_BACK, UP_LEFT_FRONT, UP_LEFT_BACK, DOWN_RIGHT_FRONT, DOWN_RIGHT_BACK, DOWN_LEFT_FRONT, DOWN_LEFT_BACK};

class OctreeNode
{
 public:
  OctreeNode();
  OctreeNode(node_types_t node_type, int depth, OctreeNode *parent, double halfDistance, Point3d center);
  ~OctreeNode();
  
  std::vector<Point3d> GetObjects(Point3d p);
  bool Insert(Point3d o);
  void PrintOctreeNode();
 private: 
  Point3d _center;
  double _halfDistance; 
  int _depth;
  node_types_t _node_type;
  OctreeNode *_children;
  OctreeNode *_parent;
  std::vector<Point3d> _objects;

  bool Intersects(Point3d p);
  bool SubdivideNode();
  bool HasChildren();
  bool HasObjects();
  bool IsIn(Point3d p);
};

#endif
