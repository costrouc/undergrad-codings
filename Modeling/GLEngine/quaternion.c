#include "quaternion.h"

void setQuat_f(Quaternion Q, float real, float i, float j, float k) {
	Q.w = real;	Q.V.x = i; Q.V.y = j; Q.V.z = k;
}

void setQuat_theta(Quaternion Q, float theta_x, float theta_y, float theta_z) {
	float cos_z_2 = cosf(0.5f*theta_z);
	float cos_y_2 = cosf(0.5f*theta_y);
	float cos_x_2 = cosf(0.5f*theta_x);

	float sin_z_2 = sinf(0.5f*theta_z);
	float sin_y_2 = sinf(0.5f*theta_y);
	float sin_x_2 = sinf(0.5f*theta_x);

	// and now compute quaternion
	Q.w   = cos_z_2*cos_y_2*cos_x_2 + sin_z_2*sin_y_2*sin_x_2;
	Q.V.x = cos_z_2*cos_y_2*sin_x_2 - sin_z_2*sin_y_2*cos_x_2;
	Q.V.y = cos_z_2*sin_y_2*cos_x_2 + sin_z_2*cos_y_2*sin_x_2;
	Q.V.z = sin_z_2*cos_y_2*cos_x_2 - cos_z_2*sin_y_2*sin_x_2;
}

Quaternion addQuat(Quaternion Q, Quaternion P) {
	Quaternion sumQ;
	setQuat_f( sumQ, Q.w+P.w, Q.V.x+P.V.x, Q.V.y+P.V.y, Q.V.z+P.V.z );
	return( sumQ );
}

Quaternion subQuat(Quaternion Q, Quaternion P) {
	Quaternion diffQ;
	setQuat_f( diffQ, Q.w-P.w, Q.V.x-P.V.x, Q.V.y-P.V.y, Q.V.z-P.V.z );
	return( diffQ );
}

Quaternion multQuat(Quaternion Q, Quaternion P) {
	Quaternion productQ;
	setQuat_f( productQ,Q.w*P.w - dotV3f( Q.V, P.V ), 
	 					Q.V.y*P.V.z - Q.V.z*P.V.y + Q.w*P.V.x + Q.V.x*P.w,
						Q.V.z*P.V.x - Q.V.x*P.V.z + Q.w*P.V.y + Q.V.y*P.w,
						Q.V.x*P.V.y - Q.V.y*P.V.x + Q.w*P.V.z + Q.V.z*P.w );
	return( productQ );
}

Quaternion normQuat(Quaternion Q) {
	Quaternion normedQ;
	float len = sqrt( Q.w*Q.w + dotV3f(Q.V,Q.V) );
	setQuat_f( normedQ, Q.w/len, Q.V.x/len, Q.V.y/len, Q.V.z/len );
	return normedQ;
}