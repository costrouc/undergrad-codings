#include "ray.h"

/** newSphere **/
object_t *newSphere(char *objtype, int code) {
   object_t *obj;
   sphere_t *sphere;
   vector_t zero = {0, 0, 0};

   /* create new sphere structure */
   obj = newShape(objtype, code);
   sphere = malloc(sizeof(sphere_t));
   assert(sphere != NULL);
   ((shape_t *)(obj->derivedObject))->derivedShape = sphere;

   /* Set default values */
   sphere->center = zero;
   sphere->radius = 0;

   return(obj);
}

/** loadSphere **/
void loadSphere(FILE *inFP, object_t *obj, char *token) {
   char *attributes[] = {"center", "radius", NULL};
   int ndx;
   int expected = 0;
   int code     = 0;
   shape_t *shape = obj->derivedObject;
   sphere_t *sphere = shape->derivedShape;

   ndx = getindex(token, attributes);

   switch (ndx) {
      case 0:
         /** center **/
         expected = 3;
         code = tRead(inFP, &sphere->center);
         break;

      case 1:
         /** radius **/
         expected = 1;
         code = fscanf(inFP, "%lf", &(sphere->radius));
         break;

      default:
         expected = 0;
         loadShape(inFP, obj, token);
   }

   /* Test for bad attribute value(s) read */
   if (code != expected) {
      fprintf(stderr,
              "ERROR[loadSphere]: Object %s\n", obj->name);
      fprintf(stderr,
              "     Error reading attributes associated with token %s\n",
              token);
      /* Terminate program */
      exit(1);
   }

};

/** completeShere **/
void completeSphere(world_t *world, object_t *obj) {
     completeShape(world, obj);
}

/** dumpSphere **/
void dumpSphere(FILE *out, object_t *obj) {
   shape_t *shape = obj->derivedObject;
   sphere_t *sphere = shape->derivedShape;

   dumpShape(out, obj);

   fprintf(out, "   center:     %7.3lf %7.3lf %7.3lf\n",
            sphere->center.x,
            sphere->center.y,
            sphere->center.z);
   fprintf(out, "   radius:     %7.3lf\n",
            sphere->radius);
}

/** sphere_hits **/
hitdata_t hitSphere(object_t *obj, point_t base, vector_t dir) {
   double t;              /* Distance                          */
   point_t vPrime;        /* Remapped view point               */
   double a, b, c;        /* Intermediate results              */
   double disc;           /* Discriminant                      */
   vector_t v2;           /* Vector from viewpoint to hitpoint */
   vector_t n2;           /* Normal to hit point               */
   shape_t *shape;        /* Pointer to shape specific data    */
   sphere_t *spherePtr;   /* Pointer to sphere data            */

   hitdata_t hitdata;

   hitdata.hit = 0;

   shape = obj->derivedObject;
   spherePtr = shape->derivedShape;

   /* Remap the viewpoint */
   vPrime = tDiff(base, spherePtr->center);

   /* Set components of quadractic equation */
   a = tDot(dir, dir);
   b = 2*tDot(vPrime, dir);
   c = tDot(vPrime, vPrime) -
             spherePtr->radius * spherePtr->radius;

   /* Compute discriminant */
   disc = b*b - 4*a*c;
  
   if (disc <= 0) {
      /* Miss */
      return(hitdata);
   }

   /* This may be a hit  -- compute the distance */
   t = (-b - sqrt(disc))/(2*a);
   hitdata.distance = t;
  
   /* Set hit point for sphere and determine the normal */
   v2 = tScale(t, dir);
   hitdata.hitpoint = tSum(base, v2);
   if (hitdata.hitpoint.z > 0) {
      /* Hit is in front of screen  -- we don't display that */
      return(hitdata);
   }

   /* Compute the normal */
   n2 = tDiff(hitdata.hitpoint, spherePtr->center);
   hitdata.normal = tUnit(n2);

   /* And return the distance to the hit point */
   hitdata.hit = 1;
   hitdata.obj = obj;
   return(hitdata);
} /* End sphere_hits */
