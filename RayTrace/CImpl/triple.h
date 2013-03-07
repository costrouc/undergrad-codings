/* Prototype statements for triple library (structure version)
 *
 * CPSC 210, Spring 2011
 *
 */

#ifndef TRIPLE_H
#define TRIPLE_H

#include <stdio.h>


/** Data type definitions **/
typedef struct triple {
   double x;
   double y;
   double z;
} triple_t;


/** Triple Library Prototype Statement **/
/* Triple Length */
double tLength(triple_t t1);

/* Scale Triple
 *
 * Returns a triple containing (t * factor)
 */
triple_t tScale(
       double factor,              /* Scaling factor */
       triple_t t                  /* Input triple   */
);

/* Unitize triple
 *
 * Returns a triple containing (t / ||t||)
 */
triple_t tUnit(
       triple_t t                 /* Input triple            */
);

/* Dot Product */
double tDot(
       triple_t t1,               /* Input triple 1 */
       triple_t t2                /* Input triple 2 */
);

/* Cross Product */
triple_t tCross(
       triple_t t1,               /* Input triple 1 */
       triple_t t2                /* Input triple 2 */
);

/* tProject */
triple_t tProject(
	triple_t normal,	/* normal vector */
	triple_t v		/* regular vector */
);

/* Triple Difference
 *
 * Returns a triple containing (t1 - t2)
 */
triple_t tDiff(
       triple_t t1,               /* First triple   */
       triple_t t2                /* Second triple  */
);
     
/* Triple Addition
 * Returns a  triple containing (t1 + t2)
 */
triple_t tSum(
       triple_t t1,               /* First triple   */
       triple_t t2                /* Second triple  */
);

/* Triple print
 * Print a label and the x, y, z values of a triple 
*/
void tPrint(
      FILE *outfile,             /* Output stream */
      char *label,               /* Message lable */
      triple_t t                 /* triple        */
);

/* Triple read
 * Read values from the specified input stream into the triple
*/
int tRead(
     FILE *infile,               /* input stream      */
     triple_t *t                 /* destination triple */
);

#endif
