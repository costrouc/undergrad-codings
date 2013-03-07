#ifndef OBJECT_H
#define OBJECT_H

/** object_t -- base object type **/
typedef struct object_type {
   char    *name;       /* left_wall, center_sphere */
   char    *type;       /* plane, sphere, ...       */
   int     code;        /* Code for object type     */

   void    *derivedObject; /* Pointer to object derived from on object_t */
   void	   (*load)();
   void	   (*complete)();
   void	   (*dump)();
} object_t;

/** object_t prototypes **/
object_t   *newObject(char *objtype, int objcode);
void loadObject(FILE *inFP, object_t *obj, char *token);
void dumpObject(FILE *out, object_t *obj);

#endif
