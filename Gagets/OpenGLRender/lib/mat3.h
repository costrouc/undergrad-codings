#ifndef MAT3_H
#define MAT3_H

#include "glmathprototypes.h"

mat3_t mat3_create(mat3_t mat);
mat3_t mat3_set(mat3_t mat, mat3_t dest);
mat3_t mat3_identity(mat3_t dest);
mat3_t mat3_transpose(mat3_t mat, mat3_t dest);
mat3_t mat3_inverse(mat3_t mat, mat3_t dest);
vec3_t mat3_multiplyVec3(mat3_t mat, vec3_t vec, vec3_t dest);
mat4_t mat3_toMat4(mat3_t mat, mat4_t dest);
void mat3_str(mat3_t mat, char *buffer);

#endif
