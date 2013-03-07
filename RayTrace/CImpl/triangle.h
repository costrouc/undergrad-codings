/* triangle.h */
#ifndef TRIANGLE_H
#define TRIANGLE_H

typedef struct triangle_type{
	vector_t orient[2];
	double length[2];
}triangle_t;
 
/** polygon function prototype statements **/
object_t *newTriangle(char *objtype, int code);
void loadTriangle(FILE *inFP, object_t *obj, char *token);
void completeTriangle(world_t *world, object_t *obj);
void dumpTriangle(FILE *outFP, object_t *obj);
hitdata_t hitTriangle(object_t *obj, point_t base, vector_t dir);

#endif
