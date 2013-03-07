#ifndef VECTOR3F_H
#define VECTOR3F_H

#include <math.h>

typedef struct {
	float x,y,z; // components
} Vector3f;

void setVector3f(Vector3f V, float a, float b, float c);

Vector3f addV3f( Vector3f V, Vector3f U );

Vector3f subV3f( Vector3f V, Vector3f U );

float dotV3f( Vector3f V, Vector3f U );

#endif