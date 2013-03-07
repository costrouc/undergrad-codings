#ifndef IMAGE_H
#define IMAGE_H

/** hit_t -- data related to an object "hit" **/
typedef struct hitdata {
   int      hit;         /* Boolean indicating if hit or not */
   object_t *obj;        /* Pointer to hit object            */
   point_t   hitpoint;   /* Hit point coordiates             */
   vector_t  normal;     /* Normal at hit point              */
   double    distance;   /* Distance to hitpoint             */
} hitdata_t;

/** Prototype statements for modules image.c and raytrace.c **/
image_t *newImage(int columns, int rows, int brightness);
image_t *makeImage(world_t *world);
pixel_t makePixel(world_t *world, int colndx, int rowndx);
double randPix(int coord);
vector_t genRay(world_t *world, int column, int row);
intensity_t rayTrace(world_t *world, point_t base, vector_t unitDir,
                 double total_dist, object_t *self);
hitdata_t findClosest(list_t *sceneList, point_t base, vector_t unitDir,
                  object_t *self);

#endif
