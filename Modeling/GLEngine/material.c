#include "material.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

void setOGLMaterial(GLenum face, Material *m)
{
  glMaterialfv(face,GL_AMBIENT,m->_ambient);
  glMaterialfv(face,GL_DIFFUSE,m->_diffuse);
  glMaterialfv(face,GL_SPECULAR,m->_specular);
  glMaterialfv(face,GL_SHININESS,m->_shininess);
}
