/* plane.c */

#include "ray.h"

/** newPlane **/
object_t *newPlane(char *objtype, int code) {
   object_t *obj;
   plane_t *plane;
   vector_t zero = {0, 0, 0};

   /* create new plane structure */
   obj = newShape(objtype, code);
   plane = malloc(sizeof(plane_t));
   assert(plane != NULL);
   ((shape_t *)(obj->derivedObject))->derivedShape = plane;

   /* Set default values */
   plane->point = zero;
   plane->normal = zero;

   return(obj);
}

/** loadPlane **/
void loadPlane(FILE *inFP, object_t *obj, char *token) {
   char *attributes[] = {"point", "normal", NULL};
   int ndx;
   int expected = 0;
   int code = 0;
   shape_t *shape = obj->derivedObject;
   plane_t *plane = shape->derivedShape;

   ndx = getindex(token, attributes);

   switch (ndx) {
      case 0:
         /* point */
         expected = 3;
         code = tRead(inFP, &plane->point);
         break;

      case 1:
         /* normal */
         expected = 3;
         code = tRead(inFP, &plane->normal);
         break;

      default:
         loadShape(inFP, obj, token);
   }

   /* Test for bad attribute value(s) read */
   if (code != expected) {
      fprintf(stderr,
              "ERROR[loadPlane]: Object %s\n", obj->name);
      fprintf(stderr,
              "     Error reading attributes associated with token %s\n",
              token);
      /* Terminate program */
      exit(1);
   }
}

/** completePlane **/
void completePlane(world_t *world, object_t *obj) {
     completeShape(world, obj);
}

/** dumpPlane **/
void dumpPlane(FILE *outFP, object_t *obj) {
   shape_t *sobj = obj->derivedObject;
   plane_t *plane = sobj->derivedShape;

   dumpShape(outFP, obj);

   fprintf(outFP, "   point:      %7.3lf %7.3lf %7.3lf\n", 
            plane->point.x,
            plane->point.y,
            plane->point.z);
   fprintf(outFP, "   normal:     %7.3lf %7.3lf %7.3lf\n", 
            plane->normal.x,
            plane->normal.y,
            plane->normal.z);
}


hitdata_t hitPlane(object_t *obj, point_t base, vector_t dir) {
    shape_t *shape =    obj->derivedObject;     /* Pointer to shape_t */
    plane_t *planePtr = shape->derivedShape; /* Pointer to plane   */
    hitdata_t hitdata;                       /* Hit data           */

    point_t  Q = planePtr->point;      /* Point data      */
    vector_t N = planePtr->normal;     /* Normal data     */
    vector_t D = dir;                  /* Direction vector*/
    point_t  V = base;                 /* Base coordinates*/
    point_t  H;                        /* Hit point                 */
    double   t;                        /* Distance                  */

    hitdata.hit = 0;

    if (tDot(N, D) == 0)
        return hitdata; // parallel
    
    t = (tDot(N, Q) - tDot(N, V))/
           tDot(N, D);
    if (t < 0)
        return hitdata; // behind me
    
    H = tScale(t, D); 
    H = tSum(H, V);
    if (H.z > 0)
        return hitdata; // between the "screen" and me
    
    hitdata.hit = 1;
    hitdata.obj = obj;
    hitdata.hitpoint = H;
    hitdata.normal = N;
    hitdata.distance = t;
    
    return hitdata;
}
