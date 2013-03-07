#include "plyobject.h"
#include "phong.h"
#include "keyinput.h"
#include "mouseinput.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include<stdio.h>

int main(int argc, char **argv)
{
  if (argc != 2)
    {
      fprintf(stderr,"Error: Input must be in form: ./program inputfile.ply\n");
      return 1;
    }
  
  initOGL(argc,argv);
  glutKeyboardFunc(getKeyboardInput);
  glutMouseFunc(getMouseInput);
  glutDisplayFunc(renderScene);
  glutVisibilityFunc(visible);
  glutMainLoop();
  return 0;
}
