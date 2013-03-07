/* window.c */

#include "ray.h"

/** newWindow **/
object_t *newWindow(char *objtype, int code) { 
   object_t *obj;
   window_t *window;

   /* create new window structure */
   obj = newObject(objtype, code);
   window = malloc(sizeof(window_t));
   assert(window != NULL);
   obj->derivedObject = window;
   window->samples=1;
   window->anaglyph=0;
   return(obj);
}

/** loadWindow **/
void loadWindow(FILE *inFP, object_t *obj, char *token) {
  /** STUBBED **/
  char *attributes[] = {"worldwidth", "worldheight","pixelwidth","viewpoint","ambient","samples","anaglyph", NULL};
  int ndx;
  int expected = 0;
  int code     = 0;
  window_t *window = obj->derivedObject;
  
  ndx = getindex(token, attributes);
  
  switch (ndx) {
  case 0:
    /** worldwidth **/
    expected = 1;
    code = fscanf(inFP,"%lf",&(window->worldWidth));
    break;
    
  case 1:
    /** worldheight **/
    expected = 1;
    code = fscanf(inFP,"%lf",&(window->worldHeight));
    break;
    
  case 2:
    /** pixelwidth **/
    expected = 1;
    code = fscanf(inFP,"%d",&(window->pixelWidth));
    break;
    
  case 3:
    /** viewpoint **/
    expected = 3;
    code = tRead(inFP,&window->viewPoint);
    break;
    
  case 4:
    /** ambient **/
    expected = 3;
    code = tRead(inFP,&window->ambient);
    break;

  case 5:
    /** samples **/
    expected = 1;
    code = fscanf(inFP,"%d",&window->samples);
    break;

  case 6:
    /** anaglyph **/
    expected = 1;
    code = fscanf(inFP,"%lf",&(window->anaglyph));
    break;
    
  default:
    expected = 0;
    loadObject(inFP, obj, token);
  }
  
  /* Test for bad attribute value(s) read */
  if (code != expected) {
    fprintf(stderr,
	    "ERROR[loadWindow]: Object %s\n", obj->name);
    fprintf(stderr,
	    "     Error reading attributes associated with token %s\n",
	    token);
    /* Terminate program */
    exit(1);
  }
}

/** completeWindow **/
void completeWindow(world_t *world, object_t *obj) {
   window_t *window = obj->derivedObject;

   world->window = obj;

   /* Compute world height from world width and pixel columns and rows */
   window->pixelHeight = window->pixelWidth * window->worldHeight/window->worldWidth;

  /* Compute the length of the pixel given for anaglyph */
  if (window->anaglyph != 0)
  window->anaglyph = window->anaglyph * window->pixelWidth/window->worldWidth;

}

/** dumpWindow **/
void dumpWindow(FILE *out, object_t *obj) {
   window_t *window = obj->derivedObject;

   dumpObject(out, obj);

   fprintf(out, "   %-20s%-6d\n", 
                                 "Pixel Width:",
                                 window->pixelWidth);
   fprintf(out, "   %-20s%-6d\n", 
                                 "Pixel Heigth:",
                                 window->pixelHeight);
   fprintf(out, "   %-20s%-6.1lf\n", 
                                 "World Width",
                                 window -> worldWidth);
   fprintf(out, "   %-20s%-6.1lf\n", 
                                 "World Heigth",
                                 window -> worldHeight);
   fprintf(out, "   %-20s%-6.1lf%-6.1lf%-6.1lf \n",
                                "viewPoint",
                                (window -> viewPoint).x,
                                (window -> viewPoint).y, 
                                (window -> viewPoint).z );
   fprintf(out, "   %-20s%-6.1lf%-6.1f%-6.1lf\n", 
                                "ambient intensity:",
                                (window -> ambient).x, 
                                (window -> ambient).y, 
                                (window -> ambient).z );
}

