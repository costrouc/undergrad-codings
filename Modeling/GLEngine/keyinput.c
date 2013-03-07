#include "keyinput.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glext.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glext.h>
#endif

#include <stdio.h>
#include <stdlib.h>

void getKeyboardInput(unsigned char key, int x, int y)
{
  //printf("You hit key: %c Number: %d\n",key, (int)key);	

		switch (key) {
		  	case 27:             /* escape */
           glDeleteBuffers(3,attrLocation);
           glDeleteTextures(1,textureLocation);
           glDeleteTextures(1,textureNormalLocation);
           glDeleteTextures(1,textureDispLocation);
           glDeleteProgram(shaderProgDiffNormDisp);
           glDeleteProgram(shaderProgPhong);
           glDeleteProgram(shaderProgNormal);
           exit(0);
           break;
	  		case ' ':
		    	if (!moving) {
		      		idle();
		      		glutPostRedisplay();
		    	}
			// object rotations
		  	case 'd':
				objRot_x+=0.25;
				break;
		  	case 'a':
				objRot_x-=0.25;
				break;
			case 'w':
				objRot_y+=0.25;
				break;
			case 's':
				objRot_y-=0.25;
				break;
			case 'q':
				objRot_z+=0.25;
				break;
			case 'e':
				objRot_z-=0.25;
				break;
			// object translations
			case 'k': // left
				trans_x-=0.25;
				break;
		  	case 'o': // up
				trans_y+=0.25;
				break;
			case ';': // right
				trans_x+=0.25;
				break;
			case 'l': // down
				trans_y-=0.25;
				break;
			case 'i': 
				trans_z+=0.25;
				break;
			case 'p':
				trans_z-=0.25;
				break;
	  	}
}

void idle(void) {
  glutPostRedisplay();
}

void visible(int state) {
  	if (state == GLUT_VISIBLE) {
    	if (moving)
      		glutIdleFunc(idle);
  	}
	else {
    	if (moving)
      		glutIdleFunc(NULL);
  	}
}
