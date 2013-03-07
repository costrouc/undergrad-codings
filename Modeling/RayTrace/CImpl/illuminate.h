#ifndef ILLUMINATE_H
#define ILLUMINATE_H

intensity_t diffuse(world_t *world,hitdata_t *hitdata);
intensity_t processLight(list_t *list,hitdata_t *hitdata,object_t *obj);

#endif
