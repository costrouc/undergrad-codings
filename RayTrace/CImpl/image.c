#define XDEBUG
#include "ray.h"

/** newImage **/
image_t *newImage(int columns, int rows, int brightness) {
   image_t *new;
   new = malloc(sizeof(image_t));
   new->rows = rows;
   new->columns = columns;
   new->brightness = brightness;
   new->image = malloc(new->rows*new->columns*sizeof(pixel_t));
   return(new);
}

/** makeImage **/
image_t *makeImage(world_t *world) {

  /** STUBBED **/
  int brightness=255;
  int rowNdx,colNdx;
  pixel_t *pixel;
  pixel_t *temp=malloc(sizeof(pixel_t));
  image_t *image;
  window_t *window=((object_t *)world->window)->derivedObject;
  
  image=newImage(window->pixelWidth,window->pixelHeight,brightness);
  pixel=image->image;
  /* The reason these for loops are read weird is because it is the way 
  that the fwrite will write the data */
  for (rowNdx=(window->pixelHeight)-1;rowNdx>=0;rowNdx--){
    for (colNdx=0;colNdx<(window->pixelWidth);colNdx++){
      *pixel=makePixel(world,colNdx,rowNdx);
      if (window->anaglyph != 0){
      	 *temp=makePixel(world,(colNdx+window->anaglyph),rowNdx);
      	 pixel->g=temp->g;
      	 pixel->b=temp->b;
      }
      pixel++;
    }
  }    
  return(image);
} /* End make_image */

/** makePixel **/
pixel_t makePixel(world_t *world, int colndx, int rowndx) {
   pixel_t pixel = {0, 0, 0};
   intensity_t intensity={0,0,0};
   int sampNdx;
   window_t *window=((object_t *)world->window)->derivedObject;

   /** STUBBED **/
   for (sampNdx=0;sampNdx<window->samples;sampNdx++){
        vector_t dir=genRay(world,colndx,rowndx);
         point_t base=((window_t *)((object_t *)world->window)->derivedObject)->viewPoint;
   
         intensity=tSum(intensity,rayTrace(world,base,dir,0.0,NULL));
   }
   intensity=tScale(1.0/window->samples,intensity);

   if (intensity.x>1.0)
     intensity.x=1.0;
   if (intensity.y>1.0)
     intensity.y=1.0;
   if (intensity.z>1.0)
     intensity.z=1.0;

   pixel.r=255.0*intensity.x;
   pixel.g=255.0*intensity.y;
   pixel.b=255.0*intensity.z;
   
   return(pixel);

} /* End make_pixel */

double randPix(int coord){
	double rNum=drand48();
	return coord+rNum-0.5;
}

