#include "ray.h"


/** genRay **/
vector_t genRay(world_t *world, int column, int row) {
   window_t *window = world->window->derivedObject;
   vector_t direction;
   

   direction.x = (randPix(column)/(double)(window->pixelWidth-1))*window->worldWidth;
   direction.y = (randPix(row)/(double)(window->pixelHeight-1))*window->worldHeight;
   direction.z = 0;

   direction = tDiff(direction, window->viewPoint);
   direction = tUnit(direction);
   return(direction);
}

/** rayTrace **/
intensity_t rayTrace(world_t *world, point_t base, vector_t unitDir,
                 double total_dist, object_t *self) {

   /** STUBBED **/
   intensity_t color = {0, 0, 0};
   hitdata_t hitdata;
   
   hitdata = findClosest(world->sceneList, base, unitDir, self);
   
   if (hitdata.hit==0)
     return(color);
   
   total_dist=total_dist+hitdata.distance;
   
   color.x=((shape_t *)hitdata.obj->derivedObject)->color.r;
   color.y=((shape_t *)hitdata.obj->derivedObject)->color.g;
   color.z=((shape_t *)hitdata.obj->derivedObject)->color.b;

   color.x=(((window_t *)world->window->derivedObject)->ambient.x)*color.x;
   color.y=(((window_t *)world->window->derivedObject)->ambient.y)*color.y;
   color.z=(((window_t *)world->window->derivedObject)->ambient.z)*color.z;
   
   color=tScale((1/255.0),color);
   color=tSum(color,diffuse(world,&hitdata));
    
   color=tScale((1/total_dist),color);

   // added reflective code
   intensity_t reflective=((shape_t *)hitdata.obj->derivedObject)->reflective;
   intensity_t intensity;
   if (tLength(reflective)!=0){
	vector_t u= tUnit(tScale(-1,unitDir));
	vector_t v;
	vector_t n= tUnit(hitdata.normal);
        
        v=tUnit(tDiff(tScale(2*tDot(u,n),n),u));
	intensity=rayTrace(world,hitdata.hitpoint,v,total_dist,hitdata.obj);
        intensity.x=intensity.x*reflective.x;
	intensity.y=intensity.y*reflective.y;
	intensity.z=intensity.z*reflective.z;
	color=tSum(color,intensity);     
   }
   
   return(color);

} /* End ray_trace */

/** findClosest **/
hitdata_t findClosest(list_t *sceneList, point_t base, vector_t unitDir,
                      object_t *self) {
  
  hitdata_t closest;
  hitdata_t hit;
  object_t *obj;
  double distance=-1;
  
  closest.hit = 0;
  listReset(sceneList);
  while ((obj=listGet(sceneList))!=NULL){
    hit=((shape_t *)obj->derivedObject)->hit(obj,base,unitDir);
    if (hit.hit==1 && hit.obj != self)
      if (hit.distance<distance || distance<0){
	closest=hit;
        distance=closest.distance;
      }
  }
  
  return(closest);
} /* End find_closest_obj */
