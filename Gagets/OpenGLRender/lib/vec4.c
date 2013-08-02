#include "vec4.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

vec4_t vec4_create(vec4_t vec) {
  vec4_t dest = calloc(sizeof(float_t), 4);

  if (vec) {
    dest[0] = vec[0];
    dest[1] = vec[1];
    dest[2] = vec[2];
    dest[3] = vec[3];
  } else {
    dest[0] = dest[1] = dest[2] = dest[3] = 0;
  }

  return dest;
}

vec4_t vec4_set(vec4_t vec, vec4_t dest) {
  dest[0] = vec[0];
  dest[1] = vec[1];
  dest[2] = vec[2];
  dest[3] = vec[3];

  return dest;
}

vec4_t vec4_add(vec4_t vec, vec4_t vec2, vec4_t dest) {
  if (!dest || vec == dest) {
    vec[0] += vec2[0];
    vec[1] += vec2[1];
    vec[2] += vec2[2];
    vec[3] += vec2[3];
    return vec;
  }

  dest[0] = vec[0] + vec2[0];
  dest[1] = vec[1] + vec2[1];
  dest[2] = vec[2] + vec2[2];
  dest[3] = vec[3] + vec2[3];
  
  return dest;
}

vec4_t vec4_subtract(vec4_t vec, vec4_t vec2, vec4_t dest) {
  if (!dest || vec == dest) {
    vec[0] -= vec2[0];
    vec[1] -= vec2[1];
    vec[2] -= vec2[2];
    vec[3] -= vec2[3];
    return vec;
  }

  dest[0] = vec[0] - vec2[0];
  dest[1] = vec[1] - vec2[1];
  dest[2] = vec[2] - vec2[2];
  dest[3] = vec[3] - vec2[3];
  return dest;
}

vec4_t vec4_multiply(vec4_t vec, vec4_t vec2, vec4_t dest) {
  if (!dest || vec == dest) {
    vec[0] *= vec2[0];
    vec[1] *= vec2[1];
    vec[2] *= vec2[2];
    vec[3] *= vec2[3];
    return vec;
  }

  dest[0] = vec[0] * vec2[0];
  dest[1] = vec[1] * vec2[1];
  dest[2] = vec[2] * vec2[2];
  dest[3] = vec[3] * vec2[3];
  return dest;
}

vec4_t vec4_negate(vec4_t vec, vec4_t dest) {
  if (!dest) { dest = vec; }

  dest[0] = -vec[0];
  dest[1] = -vec[1];
  dest[2] = -vec[2];
  dest[3] = -vec[3];
  return dest;
}

vec4_t vec4_scale(vec4_t vec, float val, vec4_t dest) {
  if (!dest || vec == dest) {
    vec[0] *= val;
    vec[1] *= val;
    vec[2] *= val;
    vec[3] *= val;
    return vec;
  }

  dest[0] = vec[0] * val;
  dest[1] = vec[1] * val;
  dest[2] = vec[2] * val;
  dest[3] = vec[3] * val;
  return dest;
}

vec4_t vec4_normalize(vec4_t vec, vec4_t dest) {
  if (!dest) { dest = vec; }

  float x = vec[0], y = vec[1], z = vec[2], w = vec[3],
    len = sqrt(x * x + y * y + z * z + w * w);

  if (!len) {
    dest[0] = 0;
    dest[1] = 0;
    dest[2] = 0;
    dest[3] = 0;
    return dest;
  } else if (len == 1) {
    dest[0] = x;
    dest[1] = y;
    dest[2] = z;
    dest[3] = w;
    return dest;
  }

  len = 1 / len;
  dest[0] = x * len;
  dest[1] = y * len;
  dest[2] = z * len;
  dest[3] = w * len;
  return dest;
}

float vec4_length(vec4_t vec) {
  float x = vec[0], y = vec[1], z = vec[2], w = vec[3];
  return sqrt(x * x + y * y + z * z + w * w);
}

float vec4_dot(vec4_t vec, vec4_t vec2) {
  return vec[0] * vec2[0] + vec[1] * vec2[1] + vec[2] * vec2[2] + vec[3] * vec2[3];
}

void vec4_str(vec4_t vec, char *buffer) {
  sprintf(buffer, "[%f, %f, %f %f]", vec[0], vec[1], vec[2], vec[3]);
}
