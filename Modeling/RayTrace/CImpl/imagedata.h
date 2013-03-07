#ifndef IMAGEDATA_H
#define IMAGEDATA_H

/** pixel_t definition **/
typedef struct pixelType {
   unsigned char r;
   unsigned char g;
   unsigned char b;
} pixel_t;

/** image_t defintion **/
typedef struct imageType {
   pixel_t *image;
   int     columns;
   int     rows;
   int     brightness;
} image_t;

#endif
