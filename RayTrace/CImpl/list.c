/* List functions */

#include "list.h"

/** listCreate -- create a new list object **/
list_t *listCreate() {
    list_t *new;
    new = malloc(sizeof(list_t));
    new->head = NULL;
    new->position = NULL;
    return new;
}


/** listAdd -- add an object to the linked list **/
void listAdd(list_t *list, void *objPtr) {
   node_t *new;
   /* Create new node */
   new = malloc(sizeof(node_t));
   new->objPtr = objPtr;

   /* Link node to head of list */
   new->next = list->head;
   list->head = new;
}

/** listGet -- return the object pointed to by the node pointed to by
               position, and then advance position to the next node in 
               the list.
**/
void *listGet(list_t *list) {
    node_t *nodePtr = list->position;
    if (nodePtr == NULL) {
       return NULL;
    }

    /* Advance to next node in list */
    list->position = list->position->next;
    return (nodePtr->objPtr);
}

/** listGetN -- return the nth from the front object **/
void *listGetN(list_t *list, int n) {
    node_t *ptr = list->head;
    if (n < 1) {
        return(NULL);
    }
    while (ptr != NULL && n>1) {
        ptr = ptr->next;
        n--;
    }
    if (ptr == NULL) {
       return(NULL);
    }
    return(ptr->objPtr);
}


/** listReset -- reset position to point to first node of list **/
void listReset(list_t *list) {
    list->position = list->head;
}

/** listMore -- test if more nodes in list **/
int listMore(list_t *list) {
   if (list->position == NULL) {
      return(0);
   }
   else {
      return(1);
   }
}
