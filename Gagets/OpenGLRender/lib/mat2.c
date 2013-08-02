#include "mat2.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

mat2_t mat2_create(mat2_t mat) {
    mat2_t dest = calloc(sizeof(float), 9);

    if (mat) {
        dest[0] = mat[0];
        dest[1] = mat[1];
        dest[2] = mat[2];
        dest[3] = mat[3];
    }

    return dest;
}

mat2_t mat2_set(mat2_t mat, mat2_t dest) {
    dest[0] = mat[0];
    dest[1] = mat[1];
    dest[2] = mat[2];
    dest[3] = mat[3];
 
    return dest;
}

mat2_t mat2_identity(mat2_t dest) {
    if (!dest) { dest = mat2_create(NULL); }
    dest[0] = 1;
    dest[1] = 0;
    dest[2] = 0;
    dest[3] = 1;

    return dest;
}

mat2_t mat2_transpose(mat2_t mat, mat2_t dest) {
    // If we are transposing ourselves we can skip a few steps but have to cache some values
    if (!dest || mat == dest) {
      float a01 = mat[1], a10 = mat[2];
          
      mat[1] = a10;
      mat[2] = a01;

      return mat;
    }

    dest[0] = mat[0];
    dest[1] = mat[2];
    dest[2] = mat[1];
    dest[3] = mat[3];

    return dest;
}

mat2_t mat2_inverse(mat2_t mat, mat2_t dest)
{
  if (!dest) {dest = mat;}

  // Cache the matrix values (makes for huge speed increases!)
  float a00 = mat[0], a01 = mat[1],
    a10 = mat[2], a11 = mat[3],
    
    d = (a00 * a11 - a01 * a10),
    invDet;

  // Calculate the determinant
  if (!d) {return NULL;}
  invDet = 1 / d;

  dest[0] = (a11)*invDet;
  dest[1] = (-a01)*invDet;
  dest[2] = (-a10)*invDet;
  dest[3] = (a00)*invDet;

  return dest;
}

vec3_t mat2_multiplyVec2(mat2_t mat, vec2_t vec, vec2_t dest)
{
  if (!dest) { dest = vec; }

  float x = vec[0], y = vec[1];

  dest[0] = mat[0] * x + mat[1] * y;
  dest[1] = mat[2] * x + mat[3] * y;
  
  return dest;
}

void mat2_str(mat2_t mat, char *buffer)
{
    sprintf(buffer, "[%f, %f, %f, %f]", mat[0], mat[1], mat[2], mat[3]);
}
