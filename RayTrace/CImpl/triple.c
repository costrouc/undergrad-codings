/* CPSC 210
 *
 * triple.c
 *
 * Vector and triple functions
 *
 */

#include <math.h>
#include <stdlib.h>
#include "triple.h"

/** tLength **/
/* Computes the length of a 3-D triple */
double tLength(triple_t t1) {
   double sum = 0.0;
   
   sum = t1.x*t1.x + t1.y*t1.y + t1.z*t1.z;
   
   return sqrt(sum);
}


/** tScale **/
/* Scales a 3-D triple by a specified factor. Returns a new triple */
triple_t  tScale(double factor, triple_t t) {
   triple_t result;
   
   result.x = t.x * factor;
   result.y = t.y * factor;
   result.z = t.z * factor;

   return result;
}


/** tUnit **/
/* Returns a unit triple in the same direction as the triple v.  */
triple_t tUnit(triple_t v) {
   triple_t result = {0, 0, 0};
   double len;

   len = tLength(v);
   if (len != 0) {
      result = tScale(1.0/len, v);
   }
   
   return result;
}


/** tDot **/
/* Computes the dot product of two 3-D triples */
double tDot(triple_t t1, triple_t t2) {
   double sum;
   
   sum = (t1.x * t2.x) + (t1.y * t2.y) + (t1.z * t2.z);
   
   return sum;
}

/** tCross **/
/** Returns the cross product of two vectors*/
triple_t tCross(triple_t t1,triple_t t2) {
  triple_t result;
  result.x=t1.y*t2.z-t1.z*t2.y;
  result.y=t1.z*t2.x-t1.x*t2.z;
  result.z=t1.x*t2.y-t1.y*t2.x;
  return result;
}

/* Returns the projection of two vectors */
triple_t tProject(triple_t normal,triple_t v){
	return (tDiff(v,tScale(tDot(normal,v),normal)));
}

/** tDiff **/
/* Returns a new triple that has the value of t1 - t2 */
triple_t tDiff(triple_t t1, triple_t t2) {
   triple_t result;

   result.x = t1.x - t2.x;
   result.y = t1.y - t2.y;
   result.z = t1.z - t2.z;
   
   return result;
}


/** tSum **/
/* Returns a new triple that has the value of t1 + t2 */
triple_t tSum(triple_t t1, triple_t t2) {
   triple_t  result;
   
   result.x = t1.x + t2.x;
   result.y = t1.y + t2.y;
   result.z = t1.z + t2.z;
   
   return result;
}

/** tPrint **/
/* Print a triple to output stream, with a label */
void tPrint(FILE *outfile, char *label, triple_t v) {
   fprintf(outfile, "%s %8.4lf %8.4lf %8.4lf\n", label, v.x, v.y, v.z);
}

/** tRead **/
/* Input values from a designated input stream into a triple */
int tRead(FILE *infile, triple_t *v) {
    return(fscanf(infile, "%lf %lf %lf",
                          &v->x,
                          &v->y,
                          &v->z));
}
