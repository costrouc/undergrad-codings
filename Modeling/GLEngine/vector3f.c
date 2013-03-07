#include "vector3f.h"

void setVector3f(Vector3f V, float a, float b, float c) {
	V.x = a; V.y = b; V.z = c;
}

Vector3f addV3f( Vector3f V, Vector3f U ) {
	Vector3f sumV;
	setVector3f( sumV, V.x+U.x, V.y+U.y, V.z+U.z );
	return sumV;
}

Vector3f subV3f( Vector3f V, Vector3f U ) {
	Vector3f diffV;
	setVector3f( diffV, V.x-U.x, V.y-U.y, V.z-U.z );
	return diffV;
}

float dotV3f( Vector3f V, Vector3f U ) {
	return( V.x*U.x + V.y+U.y + V.z*U.z );
}