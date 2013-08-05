
/* 
  Here are functions that determine accessibility of points 

*/

#include "network.h"
#include "channel.h"
#include "networkaccessibility.h"

using namespace std;
using namespace voro;


/* Sets up accessibility class, which is used to determine if a point is accessible or not
   The class need two atomic networks (one is original and another "high accuracy" one
   which has large atoms replaced by small atoms.
   The analysis is done on "inflated" atoms by r_probe_sample
   Detection of channels and accessible pockets is done using r_probe_chan */
void AccessibilityClass::setupAndFindChannels(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double r_probe_chan, double r_probe_sampl){

 r_probe = r_probe_sampl;

 highAccuracyFlag = highAccuracy;

 // Create a temporary copy of the atomic network in which each atom's radius has been increased by the probe radius
 if(highAccuracy)
   {
   atmnet->copy(&analyzedAtomNet);
   orgatmnet->copy(&orgAtomNet);
   }else{
   orgatmnet->copy(&analyzedAtomNet);
   orgatmnet->copy(&orgAtomNet);
   };
 for(unsigned int i = 0; i < orgAtomNet.atoms.size(); i++){ orgAtomNet.atoms[i].radius += r_probe; }
 for(unsigned int i = 0; i < analyzedAtomNet.atoms.size(); i++){ analyzedAtomNet.atoms[i].radius += r_probe; }

 // Calculate and store the Voronoi network for this new atomic network

 new_rad_con = (container_periodic_poly *)performVoronoiDecomp(true, &analyzedAtomNet, &vornet, advCells, false, vorcells);

// CHANNEL::findChannels(&vornet, max(0.0, r_probe_chan - r_probe), &accessInfo, &channels);

 PORE::findChannelsAndPockets(&vornet, max(0.0, r_probe_chan - r_probe), &accessInfo, &pores);

 channelMapping.resize(accessInfo.size(),-1);
 pocketMapping.resize(accessInfo.size(),-1);
 n_channels = 0; n_pockets = 0; // number of channels and pockets
 for(unsigned int i = 0; i < pores.size(); i++)
   {
   if(pores[i].dimensionality>0)
    { // Channels
     for(unsigned int j = 0; j < pores[i].nodes.size(); j++)
       {
       channelMapping[pores[i].reverseIDMappings.find(j)->second] = n_channels;
       };
     n_channels++;
    }else
     { // Pockets
     for(unsigned int j = 0; j < pores[i].nodes.size(); j++)
       {
       pocketMapping[pores[i].reverseIDMappings.find(j)->second] = n_pockets;
       };
     n_pockets++;
    };

   };

// cout << "Accessibility setup: no channels = " << n_channels << " and no pockets = " << n_pockets << "/n";

}; 


/* Remove nodes that will not be used for analysis */
void AccessibilityClass::removeOverlappedNodes(){
  // Remove all nodes from Voronoi cells that lie within the sampling sphere
  for(unsigned int i = 0; i < vorcells.size(); i++){
    vorcells[i].removeOverlappedNodes(i, &analyzedAtomNet, 0); // 0 is probe radious, 0 because atoms are already inflated
    };
};




/* checks if the provided point is in accessible volume */
pair <bool,bool> AccessibilityClass::isVPointInsideAtomAndNotAccessible(Point samplingPoint){
return isPointInsideAtomAndNotAccessible(samplingPoint, -1);
}

/* checks if the provided point is in accessible volume */
/* return true if a point is accessible */
bool AccessibilityClass::isVPointAccessible(Point samplingPoint){
pair <bool,bool> answer = isPointInsideAtomAndNotAccessible(samplingPoint, -1);

if(answer.first == false && answer.second == false) return true; else return false;

}

/* checks if thr provided point on a surface of atomID is accessible */ 
pair <bool,bool> AccessibilityClass::isSPointInsideAtomAndNotAccessible(Point samplingPoint, int atomID){
return isPointInsideAtomAndNotAccessible(samplingPoint, atomID);
}

/* Checks is provide point is inside any atom and if it is accessible */
/* atomID is an atom in the original network from where the sampling point is
   this is to check if the point is inside atom in SA calculaton */
pair <bool,bool> AccessibilityClass::isPointInsideAtomAndNotAccessible(Point samplingPoint, int atomID){

 bool inside = false, overlaps = false; // flags to state if a point is inside atom, inaccessible

 Point smplPoint; // temporary sampling point

 double newAtomX, newAtomY, newAtomZ;
 int minAtomID;
 bool foundCell = new_rad_con->find_voronoi_cell(samplingPoint[0], samplingPoint[1], samplingPoint[2], newAtomX, newAtomY, newAtomZ, minAtomID);
 if(!foundCell){
     cerr << "Error: Unable to find Voronoi cell for sampled point." << "\n"
          << "Exiting..." << "\n";
     exit(1);
    };

 tempMinDistAtomID = minAtomID; // store this information so it is not needed to be recomputed
 tempPoint = samplingPoint;

 // The routine first checks if the provided point is inside atoms


 // if in SA routine, check if the sampled point is within other atom
 if(atomID >= 0){
  if(highAccuracyFlag == false)
    {
    // old check from Thomas code (pre-high accuracy)
 
    // If point in Voronoi cell of different atom, probe-atom overlap occurs because of d^2-r^2 criterion.
    if(minAtomID != atomID)
      overlaps = true;

    }else{
    // new check if high accuracy is requested
    if(analyzedAtomNet.IDmapping[minAtomID] != atomID)
      overlaps = true;
    };
  }; // finishing check for SA


 ATOM curAtom = analyzedAtomNet.atoms[minAtomID];

 // Adjust sampling point so that it lies within the Voronoi cell of interest constructed previously
 smplPoint = (samplingPoint.add(Point(curAtom.x, curAtom.y, curAtom.z).subtract(Point(newAtomX, newAtomY, newAtomZ))));


 double minDist = calcEuclideanDistance(smplPoint[0], smplPoint[1], smplPoint[2], curAtom.x, curAtom.y, curAtom.z);
 if(minDist < curAtom.radius - 0.00000001) 
   overlaps = true;


 if(highAccuracyFlag)  // additional check for high accuracy calculations to check if a point is within the original atom
   {
   curAtom = orgAtomNet.atoms[analyzedAtomNet.IDmapping[minAtomID]];
   minDist = orgAtomNet.calcDistance(smplPoint[0], smplPoint[1], smplPoint[2], &curAtom);
   if(minDist < curAtom.radius - 0.00000001) 
     overlaps = true;
   };

 tempMinDist = minDist; // store temporary (to be used in AV within range function

 inside = overlaps;

 if(inside == true) return pair<bool,bool> (inside,overlaps); // if the point is inside an atom
                                                              // terminate and return the answer



 // If the point is outside of atoms
 // The routine then checks if the point is in accessible or inaccessible volume/surface 

 curAtom = analyzedAtomNet.atoms[minAtomID];  // making sure we look at the correct atom net (hiAcc or regular)

 // Adjust sampling point so that it lies within the Voronoi cell of interest constructed previously
 samplingPoint = (samplingPoint.add(Point(curAtom.x, curAtom.y, curAtom.z).subtract(Point(newAtomX, newAtomY, newAtomZ))));

 minDist = calcEuclideanDistance(samplingPoint[0], samplingPoint[1], samplingPoint[2], curAtom.x, curAtom.y, curAtom.z);

 // If necessary, check Voronoi nodes of cell to determine accessibility of point
 if(!overlaps){
   BASIC_VCELL vcell = vorcells[minAtomID];
   Point circCenter = Point(curAtom.x, curAtom.y, curAtom.z);
   double samplingRadius = minDist;
   Point sampleRay = Point(samplingPoint[0]-curAtom.x, samplingPoint[1]-curAtom.y, samplingPoint[2]-curAtom.z);

   // Scan the nodes in the Voronoi cell to find if line can be drawn from the node to the sampling point
   bool foundNode = false;
   needToResampleFlag = false;
   if(vcell.getNumNodes() == 0){
     cerr << "Error: Voronoi cell of sampled point does not have any nodes" << "\n"
          << "Point: " << samplingPoint[0] << " " << samplingPoint[1] << " " << samplingPoint[2] << "\n"
          << "Voronoi cell is #" << minAtomID << "\n"
          << "Please contact the source code provider." << "\n"
          << "Exiting..." << "\n";
     exit(1);
   }
   for(int k = 0; k < vcell.getNumNodes(); k++){
           Point nodePoint = vcell.getNodeCoord(k);
           bool nodeInsideSphere = (calcEuclideanDistance(nodePoint[0], nodePoint[1], nodePoint[2], circCenter[0], circCenter[1], circCenter[2]) < samplingRadius);
           bool nodeInsideOtherAtom = (vornet.nodes[vcell.getNodeID(k)].rad_stat_sphere < 0.0 );
           if(!nodeInsideSphere&&!nodeInsideOtherAtom){
             Point otherRay = samplingPoint.subtract(nodePoint);
             double dotProduct = sampleRay.dot_product(otherRay);
             if(dotProduct > 0) {
               // Angle is less than 90 degrees and so the line segment intersects twice,
               // making the path not viable
             }
             else {
               // Angle is at least 90 degrees and so the line segment interesects only once, 
               // thereby representing a viable path
               foundNode = true;
               tempNodeID = vcell.getNodeID(k);
               overlaps = !accessInfo.at(vcell.getNodeID(k));
               break;
             }
           }
    }

    // Sampling failed due to lying on Voronoi cell face and numerical inaccurarcy. 
    // Record failure, resample and notify user later
    if(!foundNode){
       resampleCount++;
       resampledInfo.push_back(pair<int, Point> (minAtomID, samplingPoint));
       needToResampleFlag = true;
       };

 };

return pair<bool,bool> (inside,overlaps);

};

