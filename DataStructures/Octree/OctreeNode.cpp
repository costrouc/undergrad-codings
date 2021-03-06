#include "OctreeNode.h"

OctreeNode::OctreeNode(){}

OctreeNode::OctreeNode(node_types_t node_type, int depth, OctreeNode *parent, double halfDistance, Point3d center)
  :_center(center)
{
  _depth = depth;
  _halfDistance = halfDistance;
  _node_type = node_type;
  _parent = parent;
  _children = NULL;
}

OctreeNode::~OctreeNode()
{
  //I think there is a bug in memory management here
}

bool OctreeNode::HasChildren()
{
  if (_children == NULL)
    {
      return false;
    }
  
  return true;
}

bool OctreeNode::HasObjects()
{
  if (_objects.size() == 0)
    {
      return false;
    }
  
  return true;
}

bool OctreeNode::IsIn(Point3d p)
{
  if (p.GetX() <= _center.GetX() + _halfDistance  && p.GetX() >= _center.GetX() - _halfDistance)
    {
      if (p.GetY() <= _center.GetY() + _halfDistance  && p.GetY() >= _center.GetY() - _halfDistance)
	{
	  if (p.GetZ() <= _center.GetZ() + _halfDistance  && p.GetZ() >= _center.GetZ() - _halfDistance)
	    {
	      return true;
	    }
	}
    }
  
  return false;
}

bool OctreeNode::Insert(Point3d p)
{
  //std::cout << "(" << p.GetX() << "," << p.GetY() << "," << p.GetZ() << ")\n";
  //std::cout << "(" << _center.GetX() << "," << _center.GetY() << "," << _center.GetZ() << ")\n";
  //std::cout << "HalfDistance: " << _halfDistance << std::endl;
  if (Intersects(p))
    {
      if (HasChildren())
	{
	  for (int i=0; i<NUM_OF_CHILDREN; ++i)
	    {
	      _children[i].Insert(p);
	    }
	}
      else
	{
	  _objects.push_back(p);
	  if (_objects.size() > MAX_OBJECTS_PER_NODE && _depth < MAX_DEPTH && _halfDistance > MIN_CONTAINER_SIZE)
	    {
	      SubdivideNode();
	    }
	}
    }
  else
    {
      //std::cerr << "Tried to insert a Point3d that did not intersect the given node\n";
      return false;
    }
  return true;
}

std::vector<Point3d> OctreeNode::GetObjects(Point3d p)
{
  if (IsIn(p))
    {
      if (HasChildren())
	{
	  for (int i=0; i<NUM_OF_CHILDREN; ++i)
	    {
	      _children[i].GetObjects(p);
	    }
	}
      else 
	{
	  return _objects;
	}
    }
  //TODO there is a bug here

  std::cerr << "Point3d is not within node\n";
  return std::vector<Point3d>();

  
}

bool OctreeNode::SubdivideNode()
{
  assert(_objects.size() > MAX_OBJECTS_PER_NODE);
  assert(_depth < MAX_DEPTH);
  assert(_halfDistance > MIN_CONTAINER_SIZE);

  _children = new OctreeNode[NUM_OF_CHILDREN];

  double halfDistance = _halfDistance/2.0;

  _children[0] = OctreeNode(UP_RIGHT_FRONT, _depth+1, this, halfDistance , Point3d(_center.GetX() + halfDistance,_center.GetY() + halfDistance,_center.GetZ() + halfDistance));

  _children[1] = OctreeNode(UP_RIGHT_BACK, _depth+1, this, halfDistance , Point3d(_center.GetX() - halfDistance,_center.GetY() + halfDistance,_center.GetZ() + halfDistance));

  _children[2] = OctreeNode(UP_LEFT_FRONT, _depth+1, this, halfDistance , Point3d(_center.GetX() + halfDistance,_center.GetY() - halfDistance,_center.GetZ() + halfDistance));

  _children[3] = OctreeNode(UP_LEFT_BACK, _depth+1, this, halfDistance , Point3d(_center.GetX() + halfDistance,_center.GetY() - halfDistance,_center.GetZ() - halfDistance));

  _children[4] = OctreeNode(DOWN_RIGHT_FRONT, _depth+1, this, halfDistance , Point3d(_center.GetX() + halfDistance,_center.GetY() + halfDistance,_center.GetZ() - halfDistance));

  _children[5] = OctreeNode(DOWN_RIGHT_BACK, _depth+1, this, halfDistance , Point3d(_center.GetX() - halfDistance,_center.GetY() + halfDistance,_center.GetZ() - halfDistance));

  _children[6] = OctreeNode(DOWN_LEFT_FRONT, _depth+1, this, halfDistance , Point3d(_center.GetX() + halfDistance,_center.GetY() - halfDistance,_center.GetZ() - halfDistance));

  _children[7] = OctreeNode(DOWN_LEFT_BACK, _depth+1, this, halfDistance , Point3d(_center.GetX() - halfDistance,_center.GetY() - halfDistance,_center.GetZ() - halfDistance));

  for (int j=0; j< (int) _objects.size(); ++j)
    {
      for (int i=0; i<NUM_OF_CHILDREN; ++i)
	{
	  _children[i].Insert(_objects.at(j));
	}
    }
  _objects.clear();
  
  return true;
}

//This is to make it work with Point3ds
bool OctreeNode::Intersects(Point3d p)  
{  
  if (p.GetX() <= _center.GetX() + _halfDistance  && p.GetX() >= _center.GetX() - _halfDistance)
    {
      if (p.GetY() <= _center.GetY() + _halfDistance  && p.GetY() >= _center.GetY() - _halfDistance)
	{
	  if (p.GetZ() <= _center.GetZ() + _halfDistance  && p.GetZ() >= _center.GetZ() - _halfDistance)
	    {
	      return true;
	    }
	}
    }
  return false;
}

void OctreeNode::PrintOctreeNode()
{
  std::cout << "=== Depth: " << _depth << " ==========";
  std::cout << "Center: (" << _center.GetX() << "," << _center.GetY() << "," << _center.GetZ() << ")\n";
  std::cout << "HalfDistance: " << _halfDistance << std::endl;
  std::cout << "Has Children? " << HasChildren() << std::endl;
  std::cout << "Has Objects? " << HasObjects() << std::endl;
  std::cout << "Has " << _objects.size() << " objects.\n" << std::flush;

  if (HasChildren())
    {
      for (int i=0; i<NUM_OF_CHILDREN; ++i)
	{
	  _children[i].PrintOctreeNode();
	}
    }
}
