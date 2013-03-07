#ifndef WINDOW_H
#define WINDOW_H

/** window_t -- data related to scene environment **/
typedef struct window_type
{
   /** Input values **/
   double worldWidth;      /* Screen width in world coordinates              */
   double worldHeight;     /* Screen height in world coordiantes             */
   point_t  viewPoint;     /* Viewpt Loc in world coords                     */
   int      pixelWidth;    /* Pixel columns                                  */
   intensity_t ambient;    /* Ambient light level                	     */
   int samples;  	   /* Number of Samples of each pixel Anti- Aliasing */
   double anaglyph;        /* This allows 3D features */

   /** Computed values **/
   int     pixelHeight;    /* Pixel rows                         */
}  window_t;

/** Window function prototypes **/
object_t *newWindow(char *objtype, int objcode);
void loadWindow(FILE *inFP, object_t *obj, char *token);
void completeWindow(world_t *world, object_t *obj);
void dumpWindow( FILE   *out, object_t *obj);
void setWindow(world_t *world, int pixelWidth);

#endif
