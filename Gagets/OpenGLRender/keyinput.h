#ifndef KEYINPUT_H
#define HEYINPUT_H

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

void getKeyboardInput(unsigned char key, int x,int y);
void idle(void);
void visible(int state);

extern GLuint *attrLocation;
extern GLuint *textureLocation;
extern GLuint *textureNormalLocation;
extern GLuint *textureDispLocation;
extern GLuint shaderProgDiffNormDisp;
extern GLuint shaderProgPhong;
extern GLuint shaderProgNormal;

extern float objRot_x, objRot_y, objRot_z;
extern float trans_x, trans_y, trans_z;
extern int moving;

#endif
