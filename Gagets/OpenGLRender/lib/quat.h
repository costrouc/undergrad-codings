#ifndef QUAT_H
#define QUAT_H

#include "glmathprototypes.h"

quat_t quat_create(quat_t quat);
quat_t quat_set(quat_t quat, quat_t dest);
quat_t quat_calculateW(quat_t quat, quat_t dest);
float quat_dot(quat_t quat, quat_t quat2);
quat_t quat_inverse(quat_t quat, quat_t dest);
quat_t quat_conjugate(quat_t quat, quat_t dest);
float quat_length(quat_t quat);
quat_t quat_normalize(quat_t quat, quat_t dest);
quat_t quat_multiply(quat_t quat, quat_t quat2, quat_t dest);
quat_t quat_multiplyVec3(quat_t quat, vec3_t vec, vec3_t dest);
mat3_t quat_toMat3(quat_t quat, mat3_t dest);
quat_t quat_toMat4(quat_t quat, mat4_t dest);
quat_t quat_slerp(quat_t quat, quat_t quat2, float slerp, quat_t dest);
void quat_str(quat_t quat, char *buffer);

#endif
