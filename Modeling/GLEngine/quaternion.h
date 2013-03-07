#ifndef QUATERNION_H
#define QUATERNION_H

#include "vector3f.h"

typedef struct {
	float w; // real part
	Vector3f V; // imaginary vector
} Quaternion;

void setQuat_f(Quaternion Q, float real, float i, float j, float k);

void setQuat_theta(Quaternion Q, float theta_x, float theta_y, float theta_z);

Quaternion addQuat(Quaternion Q, Quaternion P);

Quaternion subQuat(Quaternion Q, Quaternion P);

Quaternion multQuat(Quaternion Q, Quaternion P);

Quaternion normQuat(Quaternion Q);

#endif