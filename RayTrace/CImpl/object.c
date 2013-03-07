#include "ray.h"

object_t *newObject(char *type, int code) {
   object_t *new = (object_t *)malloc(sizeof(object_t));

   assert(new != NULL);

   new->name = NULL;
   new->type =  strdup(type);
   new->code = code;

    void (*loadAttribute[7])(FILE *inFP,object_t *obj,char *attribute);
    loadAttribute[0]=&loadWindow;
    loadAttribute[1]=&loadPlane;
    loadAttribute[2]=&loadSphere;
    loadAttribute[3]=&loadPolygon;
    loadAttribute[4]=&loadLight;
    loadAttribute[5]=&loadTriangle;
    loadAttribute[6]=&loadRectangle;
   
    void (*completeTypes[7])(world_t *world,object_t *obj);
    completeTypes[0]=&completeWindow;
    completeTypes[1]=&completePlane;
    completeTypes[2]=&completeSphere;
    completeTypes[3]=&completePolygon;
    completeTypes[4]=&completeLight;
    completeTypes[5]=&completeTriangle;
    completeTypes[6]=&completeRectangle;

    void (*dump[7])(FILE *outFP, object_t *obj);
    dump[0]=&dumpWindow;
    dump[1]=&dumpPlane;
    dump[2]=&dumpSphere;
    dump[3]=&dumpPolygon;
    dump[4]=&dumpLight;
    dump[5]=&dumpTriangle;
    dump[6]=&dumpRectangle;

    
    new->load=loadAttribute[code];
    new->complete=completeTypes[code];
    new->dump=dump[code];
       
    return(new);
}

void loadObject(FILE *inFP, object_t *obj, char *token) {
   int count;
   char inName[32];

   if (strcasecmp(token, "name") != 0) {
      fprintf(stderr, "Unknown attribute \"%s\"\n", token);
      exit(1);
   }

   /* Read the object's name */
   count = fscanf(inFP, "%31s", inName);
   assert(count==1);
   obj->name = strdup(inName);
}

void dumpObject(FILE *outFP, object_t *obj) {
   fprintf(outFP, "\nType: %s\n", obj->type);
   if (obj->name != NULL) {
      fprintf(outFP, "   Name: %s\n", obj->name);
   }
}

