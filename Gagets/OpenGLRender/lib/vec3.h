#ifndef VEC3_H
#define VEC3_H

#include "glmathprototypes.h"

vec3_t vec3_create(vec3_t vec);
vec3_t vec3_set(vec3_t vec, vec3_t dest);
vec3_t vec3_add(vec3_t vec, vec3_t vec2, vec3_t dest);
vec3_t vec3_subtract(vec3_t vec, vec3_t vec2, vec3_t dest);
vec3_t vec3_multiply(vec3_t vec, vec3_t vec2, vec3_t dest);
vec3_t vec3_negate(vec3_t vec, vec3_t dest);
vec3_t vec3_scale(vec3_t vec, float val, vec3_t dest);
vec3_t vec3_normalize(vec3_t vec, vec3_t dest);
vec3_t vec3_cross (vec3_t vec, vec3_t vec2, vec3_t dest);
float vec3_length(vec3_t vec);
float vec3_dot(vec3_t vec, vec3_t vec2);
vec3_t vec3_direction (vec3_t vec, vec3_t vec2, vec3_t dest);
vec3_t vec3_lerp(vec3_t vec, vec3_t vec2, float lerp, vec3_t dest);
float vec3_dist(vec3_t vec, vec3_t vec2);
vec3_t vec3_unproject(vec3_t vec, mat4_t view, mat4_t proj, vec4_t viewport, vec3_t dest);
void vec3_str(vec3_t vec, char *buffer);

#endif
