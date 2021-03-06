#ifndef UTILS_H
#define UTILS_H

#define MAXLINELENGTH 1024
#define MAXTOKENLENGTH 256

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define BUFFER_OFFSET(i) ((char *)NULL + (i))

void checkGLErrors(char *label);

#endif
