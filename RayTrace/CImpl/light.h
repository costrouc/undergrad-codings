#ifndef LIGHT_H
#define LIGHT_H

/** light_t -- light objects **/
typedef struct lightType {
	intensity_t color; /* Light rgb color */
	point_t center;	/* Center of the light */
	void *derivedLight; /* Derived object data */
} light_t;


/** Light function prototypes **/
object_t *newLight(char *objtype, int objcode);
void loadLight(FILE *inFP, object_t *obj, char *token);
void completeLight(world_t *world, object_t *obj);
void dumpLight( FILE   *out, object_t *obj);

#endif
