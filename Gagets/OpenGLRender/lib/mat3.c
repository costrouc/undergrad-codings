#include "mat3.h"
#include "mat4.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

mat3_t mat3_create(mat3_t mat) {
    mat3_t dest = calloc(sizeof(float), 9);

    if (mat) {
        dest[0] = mat[0];
        dest[1] = mat[1];
        dest[2] = mat[2];
        dest[3] = mat[3];
        dest[4] = mat[4];
        dest[5] = mat[5];
        dest[6] = mat[6];
        dest[7] = mat[7];
        dest[8] = mat[8];
    }

    return dest;
}

mat3_t mat3_set(mat3_t mat, mat3_t dest) {
    dest[0] = mat[0];
    dest[1] = mat[1];
    dest[2] = mat[2];
    dest[3] = mat[3];
    dest[4] = mat[4];
    dest[5] = mat[5];
    dest[6] = mat[6];
    dest[7] = mat[7];
    dest[8] = mat[8];
    return dest;
}

mat3_t mat3_identity(mat3_t dest) {
    if (!dest) { dest = mat3_create(NULL); }
    dest[0] = 1;
    dest[1] = 0;
    dest[2] = 0;
    dest[3] = 0;
    dest[4] = 1;
    dest[5] = 0;
    dest[6] = 0;
    dest[7] = 0;
    dest[8] = 1;
    return dest;
}

mat3_t mat3_transpose(mat3_t mat, mat3_t dest) {
    // If we are transposing ourselves we can skip a few steps but have to cache some values
    if (!dest || mat == dest) {
        float a01 = mat[1], a02 = mat[2],
            a12 = mat[5];

        mat[1] = mat[3];
        mat[2] = mat[6];
        mat[3] = a01;
        mat[5] = mat[7];
        mat[6] = a02;
        mat[7] = a12;
        return mat;
    }

    dest[0] = mat[0];
    dest[1] = mat[3];
    dest[2] = mat[6];
    dest[3] = mat[1];
    dest[4] = mat[4];
    dest[5] = mat[7];
    dest[6] = mat[2];
    dest[7] = mat[5];
    dest[8] = mat[8];
    return dest;
}

mat3_t mat3_inverse(mat3_t mat, mat3_t dest)
{
  if (!dest) {dest = mat;}

  // Cache the matrix values (makes for huge speed increases!)
  float a00 = mat[0], a01 = mat[1], a02 = mat[2],
    a10 = mat[3], a11 = mat[4], a12 = mat[5],
    a20 = mat[6], a21 = mat[7], a22 = mat[8],

    b01 = a22 * a11 - a12 * a21,
    b11 = -a22 * a10 + a12 * a20,
    b21 = a21 * a10 - a11 * a20,

    d = a00 * b01 + a01 * b11 + a02 * b21,
    id;

  if (!d) { return NULL; }
  id = 1 / d;
  
  if (!dest) { dest = mat3_create(NULL); }

  dest[0] = b01 * id;
  dest[1] = (-a22 * a01 + a02 * a21) * id;
  dest[2] = (a12 * a01 - a02 * a11) * id;
  dest[3] = b11 * id;
  dest[4] = (a22 * a00 - a02 * a20) * id;
  dest[5] = (-a12 * a00 + a02 * a10) * id;
  dest[6] = b21 * id;
  dest[7] = (-a21 * a00 + a01 * a20) * id;
  dest[8] = (a11 * a00 - a01 * a10) * id;

  return dest;
}

vec3_t mat3_multiplyVec3(mat3_t mat, vec3_t vec, vec3_t dest)
{
  if (!dest) { dest = vec; }

  float x = vec[0], y = vec[1], z = vec[2];

  dest[0] = mat[0] * x + mat[1] * y + mat[2] * z;
  dest[1] = mat[3] * x + mat[4] * y + mat[5] * z;
  dest[2] = mat[6] * x + mat[7] * y + mat[9] * z;

  return dest;
}

mat4_t mat3_toMat4(mat3_t mat, mat4_t dest) {
    if (!dest) { dest = mat4_create(NULL); }

    dest[15] = 1;
    dest[14] = 0;
    dest[13] = 0;
    dest[12] = 0;

    dest[11] = 0;
    dest[10] = mat[8];
    dest[9] = mat[7];
    dest[8] = mat[6];

    dest[7] = 0;
    dest[6] = mat[5];
    dest[5] = mat[4];
    dest[4] = mat[3];

    dest[3] = 0;
    dest[2] = mat[2];
    dest[1] = mat[1];
    dest[0] = mat[0];

    return dest;
}

void mat3_str(mat3_t mat, char *buffer)
{
    sprintf(buffer, "[%f, %f, %f, %f, %f, %f, %f, %f, %f]", mat[0], mat[1], mat[2], mat[3], mat[4], mat[5], mat[6], mat[7], mat[8]);
}
