#ifndef WORLD_H
#define WORLD_H

/** world_t -- pointers to scene and window data **/
typedef struct world_type
{
   object_t *window;       /* Pointer to window object        */
   list_t   *sceneList;    /* "Shapes" list                   */
   list_t   *lightList;    /* "Lights" list                   */
}  world_t;

/** world prototype statements **/
int getindex(char *token, char *attributes[]);
world_t *newWorld();
void loadWorld(FILE *inFP, world_t *world);
void dumpWorld(FILE *out, world_t *world);


#endif
