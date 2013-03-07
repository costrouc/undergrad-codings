/* rectangle.h */
#ifndef RECTANGLE_H
#define RECTANGLE_H

typedef struct rectangle_type{
	vector_t orient;
	double dimension[2];
}rectangle_t;

/** rectangle function prototype statements **/
object_t *newRectangle(char *objtype, int code);
void loadRectangle(FILE *inFP, object_t *obj, char *token);
void completeRectangle(world_t *world, object_t *obj);
void dumpRectangle(FILE *outFP, object_t *obj);
hitdata_t hitRectangle(object_t *obj, point_t base, vector_t dir);

#endif
