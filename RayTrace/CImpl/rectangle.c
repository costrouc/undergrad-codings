/* rectangle.c */

#include "ray.h"

/** newRectangle **/
object_t *newRectangle(char *objtype, int code) {
   object_t *obj;
   edge_t *new;

   obj = newPolygon(objtype, code);

   new = malloc(4*sizeof(edge_t));
   shape_t *shape=obj->derivedObject;
   plane_t *plane=shape->derivedShape;
   polygon_t *polygon=plane->derivedPlane;
   polygon->edges=new;
   polygon->numEdges=4;
   polygon->derivedPolygon=malloc(sizeof(rectangle_t));
   return(obj);
}

/** loadRectangle **/
void loadRectangle(FILE *inFP, object_t *obj, char *token) {
   char *attributes[] = {"orient", "dimensions", NULL};
   int ndx;
   int expected = 0;
   int code = 0;
   shape_t *shape = obj->derivedObject;
   plane_t *plane = shape->derivedShape;
   polygon_t *polygon = plane->derivedPlane;
   rectangle_t *rectangle=polygon->derivedPolygon;
   ndx = getindex(token, attributes);

   switch (ndx) {
      case 0:
         /* point */
         expected = 3;
         code = tRead(inFP, &(rectangle->orient));	
         break;

      case 1:
         /* normal */
         expected = 2;
         code = fscanf(inFP,"%lf %lf",&(rectangle->dimension[0]),&(rectangle->dimension[1]));
         break;

      default:
         loadPolygon(inFP, obj, token);
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

/** completeRectangle **/
void completeRectangle(world_t *world, object_t *obj) {
  vector_t normal;
  vector_t result;
  shape_t *shape = obj->derivedObject;
  plane_t *plane = shape->derivedShape;
  polygon_t *polygon = plane->derivedPlane;
  edge_t *edges=polygon->edges;
  rectangle_t *rectangle=polygon->derivedPolygon;
  
  normal=tUnit(plane->normal);

  result=tUnit(tProject(normal,tUnit(rectangle->orient)));
  edges[1].point=tSum(plane->point,tScale(rectangle->dimension[0],result));

  result=tUnit(tCross(normal,tUnit(rectangle->orient)));
  result=tScale(rectangle->dimension[1],result);
  edges[2].point=tSum(edges[1].point,result);
  edges[3].point=tSum(plane->point,result);

  completePolygon(world, obj);
}

/** dumpRectangle **/
void dumpRectangle(FILE *outFP, object_t *obj) {
  dumpPolygon(outFP,obj);
}

/* hit Rectangle */
hitdata_t hitRectangle(object_t *obj, point_t base, vector_t dir) {
     return hitPolygon(obj,base,dir);
}
