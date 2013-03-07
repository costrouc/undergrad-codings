/* shape.c */

#include "ray.h"

/** newShape **/
object_t *newShape(char *type, int code) {
   object_t *obj;                   /* Pointer to new object_t */
   shape_t *new;                    /* Pointer to new shape_t  */
   pixel_t white = {255, 255, 255}; /* Default color           */
   intensity_t zero = {0, 0, 0};    /* Default diffuse/reflective values */

   /* create new shape structure */
   obj = newObject(type, code);
   new = malloc(sizeof(shape_t));
   assert(new != NULL);
   obj->derivedObject = new;

   /* Set default values for a shape */
   new->color = white;
   new->diffuse = zero;
   new->reflective = zero;

   /* Set the appropriate hit function for shape_t that it is associated with */
   hitdata_t (*hitFunc[7])(object_t *obj, point_t base, vector_t dir);
   hitFunc[1]=&hitPlane;
   hitFunc[2]=&hitSphere;
   hitFunc[3]=&hitPolygon;
   hitFunc[5]=&hitTriangle;
   hitFunc[6]=&hitRectangle;

   new->hit=hitFunc[code];

   return(obj);
}


/** loadShape **/
void loadShape(FILE *inFP, object_t *obj, char *token) {
   /* Attributes recognized by loadShape */
   char *attributes[] = {"color", "diffuse", "reflective", NULL};
   int attributeNdx;    /* Index from attributes array that matches token */
   int code=0;          /* Success/failure code from processing values    */
   int expected=0;      /* Expected code result from processing values    */
   shape_t *shape = obj->derivedObject; /* Pointer to object's shape_t       */

   /* Get the index of the attribute */
   attributeNdx = getindex(token, attributes);

   /* Process the attribute */
   switch (attributeNdx) {
      case 0:
         /** color  attribute **/
         expected = 1;                          /* Success code value */
         code = getColor(inFP, &shape->color);
         break;

      case 1:
         /** diffuse attribute **/
         expected = 3;
         code = tRead(inFP, &shape->diffuse);
         break;

      case 2:
         /** reflective attribute **/
         expected = 3;
         code = tRead(inFP, &shape->reflective);
         break;

      default:
         expected = 0;
         loadObject(inFP, obj, token);
   }

   /* Test for bad attribute value(s) read */
   if (code != expected) {
      fprintf(stderr, 
              "ERROR[loadShape]: Object %s\n", obj->name);
      fprintf(stderr, 
              "     Error reading attributes associated with token %s\n",
              token);
      /* Terminate program */
      exit(1);
   }
}

/** completeShape **/
void completeShape(world_t *world, object_t *obj) {
   /* Just add new shape to sceneList linked list */
   listAdd(world->sceneList, obj);
}

/** dumpShape **/
void dumpShape(FILE *out, object_t *obj) {
   shape_t *shape = obj->derivedObject;

   /* First dump the object_t data */
   dumpObject(out, obj);

   /* And now the shape_t specific data */
   fprintf(out, "   color:       %6d  %6d  %6d\n", 
                                    shape->color.r, 
                                    shape->color.g,
                                    shape->color.b);
   tPrint(out,  "   diffuse:   ", shape->diffuse);
   tPrint(out,  "   reflective:", shape->reflective);
}
