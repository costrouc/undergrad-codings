/** main -- ray tracer 1 **/
#include "ray.h"

int main( int argc, char *argv[])
{
   /** STUBBED **/
   FILE *inFile;
   image_t *outImage;
   FILE *outFile;
   world_t *world;


   /* Test and fetch input parameters */
   if (argc != 3) {
      fprintf(stderr, "Usage: ./ray1 intxtfile outpmpfile\n");
      exit(1);
   }

   if ((outFile = fopen(argv[2], "w")) == NULL) {
      fprintf(stderr, "Cannot open output file\n");
      exit(1);
   }

   /* Input the source MDL */
   if ((inFile = fopen(argv[1],"r")) == NULL) {
      fprintf(stderr, "open MLD file failed\n");
      exit(1);
   }
   
   /*These functions do some rather dirty work more information can be found in the 
    functions */
   world=newWorld();
   printf("Created World\n");
   loadWorld(inFile,world);
   printf("Loaded World\n");
   dumpWorld(stdout,world);
   printf("Dumped World\n");
   outImage=makeImage(world);
   printf("Made Image\n");

   /* And output the image */
   fprintf(outFile, "P6 %d %d %d\n", 
     outImage->columns, outImage->rows, outImage->brightness);
   fwrite(outImage->image, sizeof(pixel_t), 
            outImage->columns * outImage->rows, outFile);
	
   exit(0);
}
