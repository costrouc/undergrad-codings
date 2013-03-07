/* polygon.c */

#include "ray.h"

/** newPolygon **/
object_t *newPolygon(char *objtype, int code) {
   object_t *obj;
   polygon_t *new;

   obj = newPlane(objtype, code);

   new = malloc(sizeof(polygon_t));
   shape_t *shape=obj->derivedObject;
   plane_t *plane=shape->derivedShape;
   plane->derivedPlane=new;
   return(obj);
}

/** loadPolygon **/
void loadPolygon(FILE *inFP, object_t *obj, char *token) {
   char *attributes[] = {"numedges", "polypoint", NULL};
   int ndx;
   int expected = 0;
   int code = 0;
   int edgeNdx;
   shape_t *shape = obj->derivedObject;
   plane_t *plane = shape->derivedShape;
   polygon_t *polygon = plane->derivedPlane;
   ndx = getindex(token, attributes);

   switch (ndx) {
      case 0:
         /* point */
         expected = 1;
         code = fscanf(inFP,"%d",&polygon->numEdges);
         polygon->edges=malloc((polygon->numEdges)*sizeof(edge_t));	
         break;

      case 1:
         /* normal */
         expected = 3;
         code= fscanf(inFP,"%d",&edgeNdx);
         code = tRead(inFP, &(polygon->edges[edgeNdx].point));
         break;

      default:
         loadPlane(inFP, obj, token);
   }

   /* Test for bad attribute value(s) read */
   if (code != expected) {
      fprintf(stderr,
              "ERROR[loadPolygon]: Object %s\n", obj->name);
      fprintf(stderr,
              "     Error reading attributes associated with token %s\n",
              token);
      /* Terminate program */
      exit(1);
   }
}

/** completePlane **/
void completePolygon(world_t *world, object_t *obj) {
  int count=0;
  completePlane(world, obj);
  shape_t *shape=obj->derivedObject;
  plane_t *plane=shape->derivedShape;
  polygon_t *polygon=plane->derivedPlane;
  polygon->edges[0].point=plane->point;

  while (count < (polygon->numEdges - 1)){
    polygon->edges[count].edgeDir=tUnit(tDiff(polygon->edges[count+1].point,polygon->edges[count].point));
    count++;
  }
  polygon->edges[count].edgeDir=tUnit(tDiff(polygon->edges[0].point,polygon->edges[count].point));     
}

/** dumpPlane **/
void dumpPolygon(FILE *outFP, object_t *obj) {
  shape_t *shape=obj->derivedObject;
  plane_t *plane=shape->derivedShape;
  polygon_t *polygon=plane->derivedPlane;
  edge_t *edge;
  int ndx=0;
  
  dumpPlane(outFP, obj);
  
  fprintf(outFP,"   Polygon data:\n");
  fprintf(outFP,"\tNumber of edges: %d\n",polygon->numEdges);
  while (ndx < (polygon->numEdges)){
    fprintf(outFP, "\tEdge %d: ",ndx);
    edge=&(polygon->edges[ndx]);
    
    fprintf(outFP, "Point:         %7.3lf %7.3lf %7.3lf\n", 
            edge->point.x,
            edge->point.y,
            edge->point.z);
    fprintf(outFP, "\t\tDirection:     %7.3lf %7.3lf %7.3lf\n", 
            edge->edgeDir.x,
            edge->edgeDir.y,
            edge->edgeDir.z);
   ndx++;
  }
}
hitdata_t hitPolygon(object_t *obj, point_t base, vector_t dir) {
	hitdata_t hit;  
	int ndx;
	shape_t *shape=obj->derivedObject;
        plane_t *plane=shape->derivedShape;
        polygon_t *polygon=plane->derivedPlane;
	edge_t *edge;
	vector_t v1,v2;
	vector_t b;
	

	hit=hitPlane(obj,base,dir);
	if (hit.hit==0)
		return hit;
	
	for(ndx=0;ndx<(polygon->numEdges);ndx++){
		edge=&(polygon->edges[ndx]);
		b=tUnit(tDiff(hit.hitpoint,edge->point));
		v1=tUnit(tCross(edge->edgeDir,b));
		if (tLength(v1)==0){
			hit.hit=0;
			return hit;
		}
		if (ndx != 0){
			if (tDot(v1,v2)<0){
				hit.hit=0;
				return hit;
			}
		}
		v2=v1;
	}
	return hit;
}
