#ifndef MAT2_H
#define MAT2_H

#include "glmathprototypes.h"

mat2_t mat2_create(mat2_t mat);
mat2_t mat2_set(mat2_t mat, mat2_t dest);
mat2_t mat2_identity(mat2_t dest);
mat2_t mat2_transpose(mat2_t mat, mat2_t dest);
mat2_t mat2_inverse(mat2_t mat, mat2_t dest);
vec2_t mat2_multiplyVec2(mat2_t mat, vec2_t vec, vec2_t dest);
void mat2_str(mat2_t mat, char *buffer);

#endif
