#ifndef SHAPE_H
#define SHAPE_H

#include "object.h"

/** shape_t -- Visible scene objects **/
typedef struct shape_type {
   /* Surface data */
   pixel_t color;            /* surface color                  */
   intensity_t  diffuse;     /* light reflection               */
   intensity_t  reflective;  /* ray reflection                 */

   void    *derivedShape;    /* Pointer to derived shape data */
   hitdata_t    (*hit)();

} shape_t;

/** shape prototype statements **/
object_t *newShape(char *type, int code);
void loadShape(FILE *inFP, object_t *obj, char *token);
void completeShape(world_t *world, object_t *obj);
void dumpShape(FILE *out, object_t *obj);

#endif
