#include "ray.h"

/** getColor()

    Converts the color name or hex value to the 3 byte RGB field

    Return 1 on success, 0 on failure.  If successfull sets r, g and
      b fields of pixel pointed to by the parameter pixel
**/
int getColor(FILE *infile, pixel_t *pixel) {
   int ndx;             /* Array index                          */
   unsigned int value;  /* RGB integer value                    */
   char word[20];       /* input word buffer                    */
   char *end;           /* set by strtol call                   */
   int code;            /* read code                            */
   /* Know color names and color values                         */
   char *colors[] =           {"white", "black", "red", "green", "blue", 
                               "orange", "purple", "yellow", NULL};
   unsigned int colorvals[] = {0xffffff, 0x000000, 0xff0000, 0x00ff00, 
                               0x0000ff, 0xff8000, 0x8000ff, 0xffff00, -1};

   /* First try to read the input as 3 decimal values */
   code = fscanf(infile, "%hhd %hhd %hhd", &pixel->r, &pixel->g, &pixel->b);
   if (code == 3) {
      return(1);     /* Success, 3 values read okay */
   }
   if (code != 0) {
      /* There must have been 1 or 2 decimal numbers, but not 3 */
      return(0);  
   }

   /* Well it's not decimal digits, is the input hex or a known color name? */
   if (fscanf(infile, "%19s", word) != 1) {
      return(0);  /* probably at end of file */
   }

   /* first check for known color names */
   if (word[0] != '#') {
      /* Search the colors array for the name */
      for (ndx=0; colors[ndx] != NULL &&
                  strcasecmp(word, colors[ndx]) != 0 ; ndx++) 
                  /* Just searching */;
      /* Note that ndx is left set to the index of the color in the
         colors array */
      if (colors[ndx] == NULL) {
         /* Color name not found */
         return(0);
      }
      /* Remember the equivalent color RGB value */
      value = colorvals[ndx];
   }
   else {
      /* Hex string of digits -- convert to internal integer*/
      value = strtol(word+1, &end, 16);
      if (*end != '\0') {
         /* Oops -- bad hex string */
         return(0);
      }
   }

   /* Finally partition the value into its three 8-bit RBG components */
   pixel->r  =  value >> 16;
   pixel->g = (value >> 8) & 0xff;
   pixel->b = value & 0xff;
   return(1);
}
