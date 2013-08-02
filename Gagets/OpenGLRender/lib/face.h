#ifndef FACE_H
#define FACE_H

#define NUMVERTPERFACE 3

/*
 *  Clipping is to be done on all of the incoming polgons such
 *  that they are broken down into individual trianges
 */

typedef struct 
{
  int _indices[3];
}face;

#endif
