#include "light.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

void setOGLLight(GLenum light,Light *l)
{
  glLightfv(light,GL_POSITION,l->_position);
  glLightfv(light,GL_AMBIENT,l->_ambient);
  glLightfv(light,GL_DIFFUSE,l->_diffuse);
  glLightfv(light,GL_SPECULAR,l->_specular);
}
