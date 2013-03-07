/** list.h -- lab 4 **/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/** List node  **/
typedef struct node_type
{
   void *objPtr;             /* Pointer to associated object  */
   struct node_type *next;   /* Pointer to next node          */
} node_t;

typedef struct list_type
{
   node_t *head;             /* Pointer to front of list      */
   node_t *position;         /* Current list position         */      
} list_t;

/** Function prototypes **/
list_t *listCreate();            /* Create and initialize list object */
void listAdd(list_t *list, void *objPtr); /* Add object to list       */
void listReset();                /* Reset position to head of list    */
void *listGet(list_t *list);     /* Get object from list              */
void *listGetN(list_t *list, int n);  /* Get nth object from list     */
int  listMore(list_t *list);     /* Test for more objects in list     */
