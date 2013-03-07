/* this is to calculate the addition to the ray */
#include "ray.h"

intensity_t diffuse(world_t *world,hitdata_t *hitdata){
  intensity_t tIlm={0,0,0};
  object_t *obj;
  
  listReset(world->lightList);
  while ((obj=listGet(world->lightList)) != NULL) {
    tIlm=tSum(tIlm,processLight(world->sceneList,hitdata,obj));
  }
  
  return tIlm;
}

intensity_t processLight(list_t *sceneList,hitdata_t *hitdata,object_t *obj){
  
  vector_t v1={0,0,0};
  vector_t ray;
  hitdata_t hit;
  vector_t result;
  vector_t dif=((shape_t *)hitdata->obj->derivedObject)->diffuse;
  light_t *light=obj->derivedObject;
  vector_t color=light->color;
  double ctheta;
  
  ray=tDiff(light->center,hitdata->hitpoint);
  
  hit=findClosest(sceneList,light->center,tUnit(tScale(-1,ray)),NULL);
  if (tDot(ray,hitdata->normal)<=0 || hit.obj !=(hitdata->obj)){
    return(v1);
  }
  
    ctheta=tDot(tUnit(hitdata->normal),tUnit(ray));
    result.x=(dif.x*color.x*ctheta)/hit.distance;
    result.y=(dif.y*color.y*ctheta)/hit.distance;
    result.z=(dif.z*color.z*ctheta)/hit.distance;
    return(result);
}      
