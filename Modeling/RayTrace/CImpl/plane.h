/* plane.h */
#ifndef PLANE_H
#define PLANE_H

/** plane_t --  infinite plane data **/
typedef struct plane_type
{
   vector_t normal;       /* Perpendicular to surface          */
   point_t  point;        /* Any point on surface              */
   void     *derivedPlane;   /* Pointer to derived object data */
}  plane_t;

/** plane function prototype statements **/
object_t *newPlane(char *objtype, int code);
void loadPlane(FILE *inFP, object_t *obj, char *token);
void completePlane(world_t *world, object_t *obj);
void dumpPlane(FILE *outFP, object_t *obj);
hitdata_t hitPlane(object_t *obj, point_t base, vector_t dir);

#endif
