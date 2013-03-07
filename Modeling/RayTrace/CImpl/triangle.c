/* triangle.c */

#include "ray.h"

/** newTriangle **/
object_t *newTriangle(char *objtype, int code) {
   object_t *obj;
   edge_t *new;

   obj = newPolygon(objtype, code);

   new = malloc(3*sizeof(edge_t));
   shape_t *shape=obj->derivedObject;
   plane_t *plane=shape->derivedShape;
   polygon_t *polygon=plane->derivedPlane;
   polygon->edges=new;
   polygon->numEdges=3;
   polygon->derivedPolygon=malloc(sizeof(triangle_t));
   return(obj);
}

/** loadTriangle **/
void loadTriangle(FILE *inFP, object_t *obj, char *token) {
   char *attributes[] = {"orient1", "orient2","length1","length2", NULL};
   int ndx;
   int expected = 0;
   int code = 0;
   shape_t *shape = obj->derivedObject;
   plane_t *plane = shape->derivedShape;
   polygon_t *polygon = plane->derivedPlane;
   triangle_t *triangle=polygon->derivedPolygon;
   ndx = getindex(token, attributes);

   switch (ndx) {
      case 0:
         /* point */
         expected = 3;
         code = tRead(inFP, &(triangle->orient[0]));	
         break;

      case 1:
         /* normal */
         expected = 3;
         code = tRead(inFP, &(triangle->orient[1]));
         break;

      case 2:
           expected = 1;
           code = fscanf(inFP,"%lf",&(triangle->length[0]));
         break;

      case 3:
           expected = 1;
           code = fscanf(inFP,"%lf",&(triangle->length[1]));
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

/** completePlane **/
void completeTriangle(world_t *world, object_t *obj) {
  vector_t normal;
  vector_t result;
  shape_t *shape = obj->derivedObject;
  plane_t *plane = shape->derivedShape;
  polygon_t *polygon = plane->derivedPlane;
  triangle_t *triangle=polygon->derivedPolygon;
  
  //I am going to compute the normal diferently b/c I don't like
  //the way we were told to compute it
  normal=tUnit(tCross(triangle->orient[0],triangle->orient[1]));
  plane->normal=normal;

  result=tUnit(tProject(normal,tUnit(triangle->orient[0])));
  polygon->edges[1].point=tSum(plane->point,tScale(triangle->length[0],result));

  result=tUnit(tProject(normal,tUnit(triangle->orient[1])));
  polygon->edges[2].point=tSum(plane->point,tScale(triangle->length[1],result));

  completePolygon(world, obj);
}

/** dumpTriangle **/
void dumpTriangle(FILE *outFP, object_t *obj) {
  dumpPolygon(outFP,obj);
}


/* hit Triangle */
hitdata_t hitTriangle(object_t *obj, point_t base, vector_t dir) {
     return hitPolygon(obj,base,dir);
}
