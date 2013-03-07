#include "ray.h"

int getindex(char *token, char *table[]) {
    int ndx = 0;
    while ((table[ndx] != NULL) && (strcasecmp(token, table[ndx]) != 0)) {
       ndx++;
    }
    return(ndx);
}

world_t *newWorld() {
   world_t *world;

   world = (world_t *)malloc(sizeof(world_t));
   world->window = NULL;
   world->sceneList = listCreate();
   world->lightList = listCreate();
   return(world);
}

void loadWorld(FILE *inFP, world_t *world) {
   object_t *obj;
   int  code;
   char *objtype=malloc(256*sizeof(char));
   char *objtypes[]={"window","plane","sphere","polygon","light","triangle","rectangle",NULL};
   char *attribute=malloc(sizeof(char)*256);
 
    object_t *(*newTypes[7])(char *objtype,int code);
    newTypes[0]=&newWindow;
    newTypes[1]=&newPlane;
    newTypes[2]=&newSphere;
    newTypes[3]=&newPolygon;
    newTypes[4]=&newLight;
    newTypes[5]=&newTriangle;
    newTypes[6]=&newRectangle;
    
    
    while(fscanf(inFP,"%s",objtype)!=EOF){
          code=getindex(objtype,objtypes);
          if (code==7){
              fprintf(stderr,"Error: Improper object type");
          exit(1);
          }

          obj=(*newTypes[code])(objtype,code);
    
          while(fscanf(inFP,"%s",attribute)){
		if (*attribute==';')
			break;               
		obj->load(inFP,obj,attribute);
          }
          obj->complete(world,obj);
    }
}

void dumpWorld(FILE *outFP, world_t *world) {

   object_t *obj;
   
   world->window->dump(outFP, world->window);

   fprintf(outFP, "\nSCENE OBJECTS:\n");
   listReset(world->sceneList);
   while ((obj=listGet(world->sceneList)) != NULL) {
      obj->dump(outFP,obj);
   }
   fprintf(outFP, "\nLIGHTS:\n");
   listReset(world->lightList);
   while ((obj=listGet(world->lightList)) != NULL) {
      obj->dump(outFP,obj);
   }
}
