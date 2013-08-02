#ifndef VEC2_H
#define VEC2_H

#include "glmathprototypes.h"

vec2_t vec2_create(vec2_t vec);
vec2_t vec2_set(vec2_t vec, vec2_t dest);
vec2_t vec2_add(vec2_t vec, vec2_t vec2, vec2_t dest);
vec2_t vec2_subtract(vec2_t vec, vec2_t vec2, vec2_t dest);
vec2_t vec2_multiply(vec2_t vec, vec2_t vec2, vec2_t dest);
vec2_t vec2_negate(vec2_t vec, vec2_t dest);
vec2_t vec2_scale(vec2_t vec, float val, vec2_t dest);
vec2_t vec2_normalize(vec2_t vec, vec2_t dest);
float vec2_length(vec2_t vec);
float vec2_dot(vec2_t vec, vec2_t vec2);
void vec2_str(vec2_t vec, char *buffer);

#endif
