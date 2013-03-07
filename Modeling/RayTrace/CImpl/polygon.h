/* polygon.h */
#ifndef POLYGON_H
#define POLYGON_H

/*list of all the points*/
typedef struct edge_type {
    point_t point;
    vector_t edgeDir;
} edge_t;


/** polygon_t -- polygon data **/
typedef struct polygon_type {
   int numEdges;
   edge_t *edges;
   void *derivedPolygon;
}  polygon_t;

/** polygon function prototype statements **/
object_t *newPolygon(char *objtype, int code);
void loadPolygon(FILE *inFP, object_t *obj, char *token);
void completePolygon(world_t *world, object_t *obj);
void dumpPolygon(FILE *outFP, object_t *obj);
hitdata_t hitPolygon(object_t *obj, point_t base, vector_t dir);

#endif
