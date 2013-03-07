#ifndef SPHERE_H
#define SPHERE_H

/** sphere_t -- sphere data **/
typedef struct sphere_type
{
   point_t    center;     /* Location of the center          */
   double     radius;     /* distance from center to surface */
   void       *derivedSphere; /* Pointer to derived object data  */
}  sphere_t;

/** sphere prototype statements **/
object_t *newSphere(char *objtype, int code);
void loadSphere(FILE *inFP, object_t *obj, char *token);
void completeSphere(world_t *world, object_t *obj);
void dumpSphere(FILE *outFP, object_t *obj);
hitdata_t  hitSphere(object_t *obj, point_t base, vector_t dir);

#endif
