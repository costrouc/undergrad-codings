#ifndef VEC4_H
#define VEC4_H

#include "glmathprototypes.h"

vec4_t vec4_create(vec4_t vec);
vec4_t vec4_set(vec4_t vec, vec4_t dest);
vec4_t vec4_add(vec4_t vec, vec4_t vec2, vec4_t dest);
vec4_t vec4_subtract(vec4_t vec, vec4_t vec2, vec4_t dest);
vec4_t vec4_multiply(vec4_t vec, vec4_t vec2, vec4_t dest);
vec4_t vec4_negate(vec4_t vec, vec4_t dest);
vec4_t vec4_scale(vec4_t vec, float val, vec4_t dest);
vec4_t vec4_normalize(vec4_t vec, vec4_t dest);
float vec4_length(vec4_t vec);
float vec4_dot(vec4_t vec, vec4_t vec2);
void vec4_str(vec3_t vec, char *buffer);

#endif
