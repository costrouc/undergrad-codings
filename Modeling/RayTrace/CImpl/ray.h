/* ray.h */
#ifndef RAY_H
#define RAY_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <memory.h>
#include <assert.h>

/** object types **/
#define WINDOW 0
#define PLANE  1
#define SPHERE 2

/** Includes of data types and function prototypes for ray tracer object.

    NOTE: be VERY careful about the order of the includes.  An object
          type (typedef) MUST be declared before it can be used, including
          by another include.  
**/
#include "triple.h"
/* For readability, equate other 3-triple types to a triple_t */
typedef triple_t vector_t;
typedef triple_t point_t;
typedef triple_t intensity_t;

#include "list.h"
#include "imagedata.h"
#include "getColor.h"
#include "object.h"
#include "world.h"
#include "image.h"
#include "shape.h"
#include "plane.h"
#include "sphere.h"
#include "light.h"
#include "window.h"
#include "polygon.h"
#include "illuminate.h"
#include "triangle.h"
#include "rectangle.h"

/** myassert macro **/
#define myassert(condition, msg, objptr) \
     if (!(condition)) { \
        fprintf(stderr, "\nERROR: %s (module %s, line %d)\n", msg, \
                __FILE__, __LINE__); \
        fprintf(stderr, "Object name %s\n", objptr->objname); \
        fprintf(stderr, "Program terminiating\n"); \
        exit(1); \
     }

#endif
