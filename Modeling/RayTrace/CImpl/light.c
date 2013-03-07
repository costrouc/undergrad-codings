/* light.c */

#include "ray.h"

/** newLight **/
object_t *newLight(char *objtype, int code) { 
   object_t *obj;
   light_t *light;

   /* create new light structure */
   obj = newObject(objtype, code);
   light = malloc(sizeof(light_t));
   assert(light != NULL);
   obj->derivedObject = light;
   return(obj);
}

/** loadWindow **/
void loadLight(FILE *inFP, object_t *obj, char *token) {
  /** STUBBED **/
  char *attributes[] = {"color","center", NULL};
  int ndx;
  int expected = 0;
  int code     = 0;
  light_t *light = obj->derivedObject;
  
  ndx = getindex(token, attributes);
  
  switch (ndx) {
  case 0:
    /** worldwidth **/
    expected = 3;
    code = tRead(inFP,&light->color);
    break;
    
  case 1:
    /** worldheight **/
    expected = 3;
    code = tRead(inFP,&light->center);
    break;
    
  default:
    expected = 0;
    loadObject(inFP, obj, token);
  }
  
  /* Test for bad attribute value(s) read */
  if (code != expected) {
    fprintf(stderr,
	    "ERROR[loadLight]: Object %s\n", obj->name);
    fprintf(stderr,
	    "     Error reading attributes associated with token %s\n",
	    token);
    /* Terminate program */
    exit(1);
  }
}

/** completeWindow **/
void completeLight(world_t *world, object_t *obj) {
   /* Just add new shape to sceneList linked list */
   listAdd(world->lightList, obj);
}

/** dumpWindow **/
void dumpLight(FILE *out, object_t *obj) {
   light_t *light = obj->derivedObject;

   dumpObject(out, obj);

   fprintf(out, "   %-20s%-6.3lf%-6.3f%-6.3lf\n", 
                                "center:",
                                (light -> center).x, 
                                (light -> center).y, 
                                (light -> center).z );
   fprintf(out, "   %-20s%-6.3lf%-6.3f%-6.3lf\n", 
                                "color:",
                                (light -> color).x, 
                                (light -> color).y, 
                                (light -> color).z );

   
}

