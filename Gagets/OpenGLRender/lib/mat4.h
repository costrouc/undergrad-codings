#ifndef MAT4_H
#define MAT4_H

#include "glmathprototypes.h"

mat4_t mat4_create(mat4_t mat);
mat4_t mat4_set(mat4_t mat, mat4_t dest);
mat4_t mat4_identity(mat4_t dest);
mat4_t mat4_transpose(mat4_t mat, mat4_t dest);
float mat4_determinant(mat4_t mat);
mat4_t mat4_toRotationMat(mat4_t mat, mat4_t dest);
mat3_t mat4_toMat3(mat4_t mat, mat3_t dest);
mat3_t mat4_toInverseMat3(mat4_t mat, mat3_t dest);
mat4_t mat4_inverse(mat4_t mat, mat4_t dest);
mat4_t mat4_multiply(mat4_t mat, mat4_t mat2, mat4_t dest);
mat4_t mat4_multiplyVec3(mat4_t mat, vec3_t vec, mat4_t dest);
mat4_t mat4_multiplyVec4(mat4_t mat, vec4_t vec, mat4_t dest);
mat4_t mat4_translate(mat4_t mat, vec3_t vec, mat4_t dest);
mat4_t mat4_scale(mat4_t mat, vec3_t vec, mat4_t dest);
mat4_t mat4_rotate(mat4_t mat, float angle, vec3_t axis, mat4_t dest);
mat4_t mat4_rotateX(mat4_t mat, float angle, mat4_t dest);
mat4_t mat4_rotateY(mat4_t mat, float angle, mat4_t dest);
mat4_t mat4_rotateZ(mat4_t mat, float angle, mat4_t dest);
mat4_t mat4_frustum(float left, float right, float bottom, float top, float near, float far, mat4_t dest);
mat4_t mat4_perspective(float fovy, float aspect, float near, float far, mat4_t dest);
mat4_t mat4_ortho(float left, float right, float bottom, float top, float near, float far, mat4_t dest);
mat4_t mat4_lookAt(vec3_t eye, vec3_t center, vec3_t up, mat4_t dest);
mat4_t mat4_fromRotationTranslation(quat_t quat, vec3_t vec, mat4_t dest);
void mat4_str(mat4_t mat, char *buffer);

#endif
