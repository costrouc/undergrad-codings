//#include "network.h"
//#include <voro++.hh>
#include <sstream>

#include "area_and_volume.h"
#include "networkinfo.h"
#include "channel.h"
#include "string_additions.h"
#include "networkaccessibility.h"
#include "general.h"

using namespace std;
using namespace voro;

/** Returns the density of the provided ATOM_NETWORK, assuming it represents a single
 *  unit cell, in g/cm^3. */
double calcDensity(ATOM_NETWORK *atmnet){
  double volume = calcDeterminant(atmnet->ucVectors); // Units of A^3
  vector<ATOM>::iterator iter = atmnet->atoms.begin();
  double massSum = 0;
  while(iter != atmnet->atoms.end()){
    massSum += iter->mass; 
    iter++;
  }
  return  massSum/(AVOGRADOS_NUMBER*volume)*1.0e24; // Units of g/cm^3
}

/* Print the coordinates contained in the provided vectors in a manner that they can be
 * displayed using ZeoVis and its VMD interface. Accessible points are colored green while
 * inaccessible points are colored red.*/
void reportPoints(ostream &output, vector<Point> axsPoints, vector<Point> inaxsPoints){
  output << "{color green}" << "\n";
  for(unsigned int i = 0; i < axsPoints.size(); i++){
    Point coords = axsPoints.at(i);
    output << "{point { " << coords[0] << " " << coords[1] << " " << coords[2] << "}}" << "\n";
  }
  
  output << "{color red}" << "\n";
  for(unsigned int i = 0; i < inaxsPoints.size(); i++){
    Point coords = inaxsPoints.at(i);
    output << "{point {" <<  coords[0] << " " << coords[1] << " " << coords[2] << "}}" << "\n";
  }
}

void reportResampledPoints(ostream &output, vector< pair<int, Point> > resampledInfo){
  output << "set num_resamples " << resampledInfo.size() << "\n";
   for(unsigned int i = 0; i < resampledInfo.size(); i++){
    Point coords = resampledInfo[i].second;
    output << "set rpoints(" << i << ") {" <<  coords[0] << " " << coords[1] << " " << coords[2] << "} " << "\n";
    output << "set rcenters(" << i << ") " << resampledInfo[i].first << "\n"; 
  }
}

/* Print the coordinates contained in the provided vectors in a manner that they can be
 *  * displayed using VisIT or analyzed using external programs using similar formating
 *  */
void reportPointsVisIT(ostream &output, vector<Point> axsPoints, vector<Point> inaxsPoints){
  for(unsigned int i = 0; i < axsPoints.size(); i++){
    Point coords = axsPoints.at(i);
    output << coords[0] << " " << coords[1] << " " << coords[2] << " 1 a" << "\n";
   }
  for(unsigned int i = 0; i < inaxsPoints.size(); i++){
    Point coords = inaxsPoints.at(i);
    output <<  coords[0] << " " << coords[1] << " " << coords[2] << " 0 n" << "\n";
  }
}

void reportPointsVisIT(ostream &output, vector<Point> axsPoints, vector<int> axsPChIDs, vector<Point> inaxsPoints, vector<int> inaxsPPIDs){
  for(unsigned int i = 0; i < axsPoints.size(); i++){
    Point coords = axsPoints.at(i);
    output << coords[0] << " " << coords[1] << " " << coords[2] << " 1 a " << axsPChIDs[i] << "\n";
   }
  for(unsigned int i = 0; i < inaxsPoints.size(); i++){
    Point coords = inaxsPoints.at(i);
    output <<  coords[0] << " " << coords[1] << " " << coords[2] << " 0 n " << inaxsPPIDs[i] << "\n";
  }
}


/* To adjust the sampling point to minimize its distance with the central atom, shift the point by the reverse of the
   atom shift vector.
 */
void adjustSamplingPoint(Point *samplingPointXYZ, Point newAtomCoordsXYZ, Point origAtomCoordsXYZ, ATOM_NETWORK *atmnet){
  Point xyzShift = origAtomCoordsXYZ.subtract(newAtomCoordsXYZ);
  Point abcShift = atmnet->xyz_to_abc(xyzShift);
  atmnet->translatePoint(samplingPointXYZ, abcShift[0], abcShift[1], abcShift[2]);
}

// Function defined in network.cc
void* performVoronoiDecomp(bool, ATOM_NETWORK *, VORONOI_NETWORK *, vector<VOR_CELL> &, bool,  vector<BASIC_VCELL> &);

/** Returns the volume accessible to a particle of the provided radius. Accessible volume is defined as any region of space
 *  in which the center of the particle of the provided radius can reach. Excludes inaccessible pockets if requested */
double calcAV(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double r_probe_chan, double r_probe, int numSamples, bool excludePockets, ostream &output, char *filename, bool visualize, bool VisITflag, bool LiverpoolFlag, bool blockingMode, double low_dist_cutoff, double high_dist_cutoff){

  //determine whether we are only sampling volume within a specific distance range of the surface
  bool within_range = false;
  if(!(low_dist_cutoff<0 || high_dist_cutoff<0)) {
    if(high_dist_cutoff<low_dist_cutoff) { //swap them
      double temp = low_dist_cutoff;
      low_dist_cutoff = high_dist_cutoff;
      high_dist_cutoff = temp;
    }
    within_range = true;
  }

  // Create an object that handles analysis of accessibility of sampled points
  AccessibilityClass accessAnalysis;
  if(highAccuracy) accessAnalysis.setupAndFindChannels(atmnet, orgatmnet, highAccuracy, r_probe_chan, r_probe);
    else accessAnalysis.setupAndFindChannels(atmnet, atmnet, highAccuracy, r_probe_chan, r_probe);

  // setting up MC sampling
  srand(randSeed);

  vector<Point> axsPoints = vector<Point> (); // stores accessible points
  vector<int> axsPointsChannelIDs; // stores the corresponding channel IDs 
  vector<Point> inaxsPoints = vector<Point> (); // inaccessible points
  vector<int> inaxsPointsPocketIDs; // stores the corresponding pocket IDs

  int count = 0;
  int count_inaxs = 0;
  int count_within_range = 0;
  int count_outside_range = 0;

  // MC statistics can be collected w.r.t to identified channels and inaccessible pockets
  vector<int> count_inChannel(accessAnalysis.n_channels,0);
  vector<int> count_inPocket(accessAnalysis.n_pockets,0);

  for(int i = 0; i < numSamples; i++){
    bool overlaps = false;
    bool inside = false;

    // Randomly sample point across the unit cell
    double aPoint = (rand()*1.0)/RAND_MAX;
    double bPoint = (rand()*1.0)/RAND_MAX;
    double cPoint = (rand()*1.0)/RAND_MAX;

    // Convert coords relative to the unit cell vectors into (x,y,z) coords
    Point samplingPoint = atmnet->abc_to_xyz(aPoint, bPoint, cPoint);

    // Calling accessibility object to determine accessibility of the point (this replaced a big chunk of code by Thomas)
    pair<bool,bool> answer = (accessAnalysis.isVPointInsideAtomAndNotAccessible(samplingPoint));
    inside = answer.first; overlaps = answer.second;
    if(accessAnalysis.needToResample() == true) i--; // the sampled point could not be analyzed in isVPointInsideAtomAndNotAccessible() function, resampling needed

    if(inside == false && excludePockets == false) overlaps = false; // if ignore inacceible pockets, treat the point as accessible (unless inside atom)
      
    // Store sampled points that did not overlap with an atom but were inaccessible for later visualization
    if(accessAnalysis.needToResample() == false && !inside && overlaps){
      count_inaxs++;
      pair <int,int> CoP = accessAnalysis.lastChannelOrPocket();
      if(CoP.first!=-1) 
        {
        cout << "Error: CoP.first!=-1 in pocket, consult source code provider\n";
        } else {
        count_inPocket[CoP.second]++;
        };
      if(!within_range) {
        Point abcCoords = Point(aPoint, bPoint, cPoint);
        Point coords = atmnet->abc_to_xyz(abcCoords);
        inaxsPoints.push_back(coords);
        inaxsPointsPocketIDs.push_back(CoP.second);
      };
    };

    // Store accessible points for later visualization
    if(accessAnalysis.needToResample() == false && !overlaps) {
      count++;     
      pair <int,int> CoP = accessAnalysis.lastChannelOrPocket();
      if(CoP.second!=-1) 
        {
        cout << "Error: CoP.second!=-1 in channel, consult source code provider\n";
        } else {
        count_inChannel[CoP.first]++;
        };
      Point abcCoords = Point(aPoint, bPoint, cPoint);
      Point coords = atmnet->abc_to_xyz(abcCoords);
      if(within_range) {
        if(accessAnalysis.lastMinDist() >= low_dist_cutoff && accessAnalysis.lastMinDist() <= high_dist_cutoff) { //NOTE: minDist is the distance between centroids, it is NOT the distance to the surface - here we are comparing this distance to the tolerance range
          count_within_range++;
          axsPoints.push_back(coords);
        } else {
          count_outside_range++;
          inaxsPoints.push_back(coords);
        }
      } else 
      {
      axsPoints.push_back(coords); 
      axsPointsChannelIDs.push_back(CoP.first);
      };
    }
  }; // ends a for loop over all sampled points
  
  // Write the necessary commands to the output stream
  // necessary to visualize the sampling
  if(visualize){
    if(VisITflag == false)
      {
      reportPoints(output, axsPoints, inaxsPoints);
      }
    else{
      if(LiverpoolFlag == false)
        {
        //report points in std. VisIt format
        reportPointsVisIT(output, axsPoints, inaxsPoints);
        }
      else{
      //report points in Liverpool format
      vector<Point> LiverpoolAxsPoints = vector<Point> ();
      vector<Point> LiverpoolInaxsPoints = vector<Point> ();
      for(unsigned i = 0; i < axsPoints.size(); i++)
         {
         LiverpoolAxsPoints.push_back(atmnet->xyz_to_abc(axsPoints[i]));
         };
      for(unsigned j = 0; j < inaxsPoints.size(); j++)
         {
         LiverpoolInaxsPoints.push_back(atmnet->xyz_to_abc(inaxsPoints[j]));
         };
      if(within_range == true) reportPointsVisIT(output, LiverpoolAxsPoints, LiverpoolInaxsPoints);
        else reportPointsVisIT(output, LiverpoolAxsPoints, axsPointsChannelIDs, LiverpoolInaxsPoints, inaxsPointsPocketIDs);
      };
      };
    //reportResampledPoints(output, resampledInfo);
  }

  // Write blocknig spheres based on current sampling points if blocking routine
  // has been requested 
  if(blockingMode)
    {
    //report points in fractonal coordinates (Liverpool format used in visualization)
    vector<Point> LiverpoolAxsPoints = vector<Point> ();
    vector<Point> LiverpoolInaxsPoints = vector<Point> ();
    for(unsigned i = 0; i < axsPoints.size(); i++)
       {
       LiverpoolAxsPoints.push_back(atmnet->xyz_to_abc(axsPoints[i]));
       };
    for(unsigned j = 0; j < inaxsPoints.size(); j++)
       {
       LiverpoolInaxsPoints.push_back(atmnet->xyz_to_abc(inaxsPoints[j]));
       };
    // call this if arg are point in factional coordiantes
    blockPockets(atmnet, output, LiverpoolAxsPoints, axsPointsChannelIDs, LiverpoolInaxsPoints, inaxsPointsPocketIDs, r_probe_chan);
    //call this if Cartesian points are to be used
//    blockPockets(atmnet, output, axsPoints, axsPointsChannelIDs, inaxsPoints, inaxsPointsPocketIDs, r_probe_chan);
    };


  // Warn user if points were resampled
  int resampleCount =  accessAnalysis.getResampleCount();
  if(resampleCount != 0){
    cerr << "\n" << "\n"
	 << "Warning: Resampled " << resampleCount << " points out of " << numSamples 
	 << " when analyzing " << atmnet->name << "\n"
	 << "\n" << "\n"; 
  }

  double volumeFraction = count*1.0/numSamples;
  double origVolume = calcDeterminant(atmnet->ucVectors);
  double av = volumeFraction * origVolume;

  // Write the results to the output stream
  if(!visualize&&!blockingMode){
    double volumeFraction_inaxs = count_inaxs*1.0/numSamples;
    double rho_crystal = calcDensity(atmnet);
    double avPerMass   = volumeFraction/rho_crystal;
    double av_inaxs = volumeFraction_inaxs * origVolume;
    double avPerMass_inaxs   = volumeFraction_inaxs/rho_crystal;
//    output << newAtomNet.name << " ";
    output << "@ " << filename << " ";
    output << "Unitcell_volume: " << origVolume << "   Density: " << rho_crystal << "   ";
    output << "AV_A^3: " << av << " "
	   << "AV_Volume_fraction: " << volumeFraction << " "
	   << "AV_cm^3/g: " << avPerMass << " "
           << "NAV_A^3: " << av_inaxs << " "
           << "NAV_Volume_fraction: " << volumeFraction_inaxs << " "
           << "NAV_cm^3/g: " << avPerMass_inaxs;
    if(within_range) {
      double volumeFraction_within_range = count_within_range*1.0/numSamples;
      double av_within_range = volumeFraction_within_range * origVolume;
      double avPerMass_within_range   = volumeFraction_within_range/rho_crystal;
      output << " range_A^3: " << av_within_range << " "
      << "range_Volume_fraction: " << volumeFraction_within_range << " "
      << "range_cm^3/g: " << avPerMass_within_range;
    }
    output << "\n";

    output << "Number_of_channels: " << count_inChannel.size()  << " Channel_volume_A^3: ";
    for(unsigned int i = 0; i < count_inChannel.size(); i++)
      {
//      output << count_inChannel[i]*100.00/count << "  ";
      output << count_inChannel[i]*(1.0/numSamples) * origVolume << "  ";
      };
    output << "\nNumber_of_pockets: " << count_inPocket.size()  << " Pocket_volume_A^3: ";
    for(unsigned int i = 0; i < count_inPocket.size(); i++)
      {
//      output << count_inPocket[i]*100.00/count_inaxs << "  ";
      output << count_inPocket[i]*(1.0/numSamples) * origVolume << "  ";
      };
    output << "\n";
  } // ends "if(!visualize){..." writing output

  accessAnalysis.deconstruct();
  return av;
}


/* backup of AV function

double calcAV(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double r_probe_chan, double r_probe, int numSamples, bool excludePockets, ostream &output, char *filename, bool visualize, bool VisITflag, bool LiverpoolFlag, double low_dist_cutoff, double high_dist_cutoff){

  //determine whether we are only sampling volume within a specific distance range of the surface
  bool within_range = false;
  if(!(low_dist_cutoff<0 || high_dist_cutoff<0)) {
    if(high_dist_cutoff<low_dist_cutoff) { //swap them
      double temp = low_dist_cutoff;
      low_dist_cutoff = high_dist_cutoff;
      high_dist_cutoff = temp;
    }
    within_range = true;
  }

  // Create a temporary copy of the atomic network in which each atom's radius has been increased by the probe radius
  ATOM_NETWORK newAtomNet;
  atmnet->copy(&newAtomNet);
  for(int i = 0; i < newAtomNet.numAtoms; i++){ newAtomNet.atoms[i].radius += r_probe; }  

  // Calculate and store the Voronoi network for this new atomic network
  VORONOI_NETWORK vornet;
  vector<BASIC_VCELL> vorcells;
  vector<VOR_CELL> advCells;

  container_periodic_poly *new_rad_con = (container_periodic_poly *)performVoronoiDecomp(true, &newAtomNet, &vornet, advCells, false, vorcells);

  vector<CHANNEL> channels = vector<CHANNEL>();
  vector<bool> accessInfo = vector<bool> ();
  CHANNEL::findChannels(&vornet, max(0.0, r_probe_chan - r_probe), &accessInfo, &channels);
  srand(randSeed);

  vector<Point> axsPoints = vector<Point> ();
  vector<Point> inaxsPoints = vector<Point> ();
  vector< pair<int, Point> > resampledInfo = vector< pair<int, Point> > ();  // List of resampled points and the id of the Voronoi cell to which they belong
  int resampleCount = 0;

  int count = 0;
  int count_inaxs = 0;
  int count_within_range = 0;
  int count_outside_range = 0;

  for(int i = 0; i < numSamples; i++){
    bool overlaps = false;

    // Randomly sample point across the unit cell
    double aPoint = (rand()*1.0)/RAND_MAX;
    double bPoint = (rand()*1.0)/RAND_MAX;
    double cPoint = (rand()*1.0)/RAND_MAX;

    // Convert coords relative to the unit cell vectors into (x,y,z) coords
    Point samplingPoint = atmnet->abc_to_xyz(aPoint, bPoint, cPoint);

    double newAtomX, newAtomY, newAtomZ;
    int minAtomID;
    bool foundCell = new_rad_con->find_voronoi_cell(samplingPoint[0], samplingPoint[1], samplingPoint[2], newAtomX, newAtomY, newAtomZ, minAtomID);
    if(!foundCell){
	cerr << "Error: Unable to find Voronoi cell for sampled point in AV calculation." << "\n"
	     << "Occurred for structure " << newAtomNet.name << "\n"
	     << "Exiting..." << "\n";
	exit(1);
      }

    ATOM curAtom = atmnet->atoms[minAtomID];

    // Adjust sampling point so that it lies within the Voronoi cell of interest constructed previously
    samplingPoint = (samplingPoint.add(Point(curAtom.x, curAtom.y, curAtom.z).subtract(Point(newAtomX, newAtomY, newAtomZ))));

    double minDist = calcEuclideanDistance(samplingPoint[0], samplingPoint[1], samplingPoint[2], curAtom.x, curAtom.y, curAtom.z);
    if(minDist < r_probe + curAtom.radius - 0.00000001)
      overlaps = true;

    bool inside = overlaps;
    
    // If necessary, check Voronoi nodes of cell to determine accessibility of point
    if(!overlaps && excludePockets){
      BASIC_VCELL vcell = vorcells[minAtomID];
      Point circCenter = Point(curAtom.x, curAtom.y, curAtom.z);
      double samplingRadius = minDist;
      Point sampleRay = Point(samplingPoint[0]-curAtom.x, samplingPoint[1]-curAtom.y, samplingPoint[2]-curAtom.z);
      
      // Scan the nodes in the Voronoi cell to find if line can be drawn from the node to the sampling point
      bool foundNode = false;
      if(vcell.getNumNodes() == 0){
	cerr << "Error: Voronoi cell of sampled point does not have any nodes" << "\n"
	     << "Point: " << samplingPoint[0] << " " << samplingPoint[1] << " " << samplingPoint[2] << "\n"
	     << "Voronoi cell is #" << minAtomID << " in structure " << newAtomNet.name << "\n"
	     << "Please contact the source code provider." << "\n"
	     << "Exiting..." << "\n";
	exit(1);
      }
      for(int k = 0; k < vcell.getNumNodes(); k++){
	      Point nodePoint = vcell.getNodeCoord(k);  
	      bool nodeInsideSphere = (calcEuclideanDistance(nodePoint[0], nodePoint[1], nodePoint[2], circCenter[0], circCenter[1], circCenter[2]) < samplingRadius);
	      if(!nodeInsideSphere){
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
	      i--;
      }
    }
      
    // Store sampled points that did not overlap with an atom but were inaccessible for later visualization
    if(!inside && overlaps){
      count_inaxs++;
      if(!within_range) {
        Point abcCoords = Point(aPoint, bPoint, cPoint);
        Point coords = atmnet->abc_to_xyz(abcCoords);
        inaxsPoints.push_back(coords);
      }
    }

    // Store accessible points for later visualization
    if(!overlaps) {
      count++;     
      Point abcCoords = Point(aPoint, bPoint, cPoint);
      Point coords = atmnet->abc_to_xyz(abcCoords);
      if(within_range) {
        if(minDist>=low_dist_cutoff && minDist<=high_dist_cutoff) { //NOTE: minDist is the distance between centroids, it is NOT the distance to the surface - here we are comparing this distance to the tolerance range
          count_within_range++;
          axsPoints.push_back(coords);
        } else {
          count_outside_range++;
          inaxsPoints.push_back(coords);
        }
      } else axsPoints.push_back(coords); 
    }
  }
  
  // Write the necessary commands to the output stream
  // necessary to visualize the sampling
  if(visualize){
    if(VisITflag == false)
      {
      reportPoints(output, axsPoints, inaxsPoints);
      }
    else{
      if(LiverpoolFlag == false)
        {
        //report points in std. VisIt format
        reportPointsVisIT(output, axsPoints, inaxsPoints);
        }
      else{
      //report points in Liverpool format
      vector<Point> LiverpoolAxsPoints = vector<Point> ();
      vector<Point> LiverpoolInaxsPoints = vector<Point> ();
      for(unsigned i = 0; i < axsPoints.size(); i++)
         {
         LiverpoolAxsPoints.push_back(atmnet->xyz_to_abc(axsPoints[i]));
         };
      for(unsigned j = 0; j < inaxsPoints.size(); j++)
         {
         LiverpoolInaxsPoints.push_back(atmnet->xyz_to_abc(inaxsPoints[j]));
         };
      reportPointsVisIT(output, LiverpoolAxsPoints, LiverpoolInaxsPoints);
      };
      };
    //reportResampledPoints(output, resampledInfo);
  }

  // Warn user if points were resampled
  if(resampleCount != 0){
    cerr << "\n" << "\n"
	 << "Warning: Resampled " << resampleCount << " points out of " << numSamples 
	 << " when analyzing " << atmnet->name << "\n"
	 << "\n" << "\n"; 
  }

  double volumeFraction = count*1.0/numSamples;
  double origVolume = calcDeterminant(atmnet->ucVectors);
  double av = volumeFraction * origVolume;

  // Write the results to the output stream
  if(!visualize){
    double volumeFraction_inaxs = count_inaxs*1.0/numSamples;
    double rho_crystal = calcDensity(atmnet);
    double avPerMass   = volumeFraction/rho_crystal;
    double av_inaxs = volumeFraction_inaxs * origVolume;
    double avPerMass_inaxs   = volumeFraction_inaxs/rho_crystal;
//    output << newAtomNet.name << " ";
    output << filename << " ";
    output << "AV_A^3: " << av << " "
	   << "AV_Volume_fraction: " << volumeFraction << " "
	   << "AV_cm^3/g: " << avPerMass << " "
           << "NAV_A^3: " << av_inaxs << " "
           << "NAV_Volume_fraction: " << volumeFraction_inaxs << " "
           << "NAV_cm^3/g: " << avPerMass_inaxs;
    if(within_range) {
      double volumeFraction_within_range = count_within_range*1.0/numSamples;
      double av_within_range = volumeFraction_within_range * origVolume;
      double avPerMass_within_range   = volumeFraction_within_range/rho_crystal;
      output << " range_A^3: " << av_within_range << " "
      << "range_Volume_fraction: " << volumeFraction_within_range << " "
      << "range_cm^3/g: " << avPerMass_within_range;
    }
    output << "\n";
  }

  delete new_rad_con;
  return av;
}

// ends backup of AV
*/

//the AV routine but instead of calculating the volume over random samples, just return a boolean array of whether or not each of a set of points is accessible
void determineAccessibility(ATOM_NETWORK *atmnet, double r_probe_chan, double r_probe, bool excludePockets, vector<bool> *isAccessible, VORONOI_NETWORK *pointSet) {
  // Create a temporary copy of the atomic network in which each atom's radius has been increased by the probe radius
  ATOM_NETWORK newAtomNet;
  atmnet->copy(&newAtomNet);
  for(int i = 0; i < newAtomNet.numAtoms; i++){ newAtomNet.atoms[i].radius += r_probe; }  

  // Calculate and store the Voronoi network for this new atomic network
  VORONOI_NETWORK vornet;
  vector<BASIC_VCELL> vorcells;
  vector<VOR_CELL> advCells;

  container_periodic_poly *new_rad_con = (container_periodic_poly *)performVoronoiDecomp(true, &newAtomNet, &vornet, advCells, false, vorcells);

  vector<CHANNEL> channels = vector<CHANNEL>();
  vector<bool> accessInfo = vector<bool> ();
  CHANNEL::findChannels(&vornet, max(0.0, r_probe_chan - r_probe), &accessInfo, &channels);
  srand(randSeed);

  vector<Point> axsPoints = vector<Point> ();
  vector<Point> inaxsPoints = vector<Point> ();
  vector< pair<int, Point> > resampledInfo = vector< pair<int, Point> > ();  // List of resampled points and the id of the Voronoi cell to which they belong
  int resampleCount = 0;

  int num_accessible = 0;

  for(unsigned int i = 0; i < pointSet->nodes.size(); i++){
    bool overlaps = false;

    Point samplingPoint;
    if(resampleCount>0) {
      double change = 0.000001;
      samplingPoint = Point(pointSet->nodes.at(i).x+change, pointSet->nodes.at(i).y+change, pointSet->nodes.at(i).z+change);
    } else samplingPoint = Point(pointSet->nodes.at(i).x, pointSet->nodes.at(i).y, pointSet->nodes.at(i).z);

    double newAtomX, newAtomY, newAtomZ;
    int minAtomID;
    bool foundCell = new_rad_con->find_voronoi_cell(samplingPoint[0], samplingPoint[1], samplingPoint[2], newAtomX, newAtomY, newAtomZ, minAtomID);
    if(!foundCell){
	    cerr << "Error: Unable to find Voronoi cell for sampled point in AV calculation." << "\n"
	    << "Occurred for structure " << newAtomNet.name << "\n"
	    << "Exiting..." << "\n";
    	exit(1);
    }

    ATOM curAtom = atmnet->atoms[minAtomID];

    // Adjust sampling point so that it lies within the Voronoi cell of interest constructed previously
    samplingPoint = (samplingPoint.add(Point(curAtom.x, curAtom.y, curAtom.z).subtract(Point(newAtomX, newAtomY, newAtomZ))));

    double minDist = calcEuclideanDistance(samplingPoint[0], samplingPoint[1], samplingPoint[2], curAtom.x, curAtom.y, curAtom.z);
    if(minDist < r_probe + curAtom.radius - 0.00000001) overlaps = true;
    
    // If necessary, check Voronoi nodes of cell to determine accessibility of point
    if(!overlaps && excludePockets){
      BASIC_VCELL vcell = vorcells[minAtomID];
      Point circCenter = Point(curAtom.x, curAtom.y, curAtom.z);
      double samplingRadius = minDist;
      Point sampleRay = Point(samplingPoint[0]-curAtom.x, samplingPoint[1]-curAtom.y, samplingPoint[2]-curAtom.z);
      
      // Scan the nodes in the Voronoi cell to find if line can be drawn from the node to the sampling point
      bool foundNode = false;
      if(vcell.getNumNodes() == 0){
	      cerr << "Error: Voronoi cell of sampled point does not have any nodes" << "\n"
	      << "Point: " << samplingPoint[0] << " " << samplingPoint[1] << " " << samplingPoint[2] << "\n"
	      << "Voronoi cell is #" << minAtomID << " in structure " << newAtomNet.name << "\n"
	      << "Please contact the source code provider." << "\n"
	      << "Exiting..." << "\n";
	      exit(1);
      }
      for(int k = 0; k < vcell.getNumNodes(); k++){
	      Point nodePoint = vcell.getNodeCoord(k);  
	      bool nodeInsideSphere = (calcEuclideanDistance(nodePoint[0], nodePoint[1], nodePoint[2], circCenter[0], circCenter[1], circCenter[2]) < samplingRadius);
	      if(!nodeInsideSphere){
	        Point otherRay = samplingPoint.subtract(nodePoint);
	        double dotProduct = sampleRay.dot_product(otherRay);
          double dotProductMax = 1e-5;
	        if(dotProduct > dotProductMax) {
	          // Angle is less than 90 degrees and so the line segment intersects twice,
      	    // making the path not viable
	        }
	        else {
	          // Angle is at least 90 degrees and so the line segment interesects only once, 
	          // thereby representing a viable path
	          foundNode = true;
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
	      i--;
      } else resampleCount = 0;
    }

    if(resampleCount==0) {
      if(!overlaps) num_accessible++;
      if(!overlaps) isAccessible->push_back(true); else isAccessible->push_back(false);
    }
  }

  // Warn user if points were resampled
  if(resampleCount != 0){
    cerr << "\n" << "\n"
    << "Warning: Resampled " << resampleCount << " points (preset coords, not random samples) out of " << pointSet->nodes.size()  
    << " when analyzing " << atmnet->name << "\n"
    << "\n" << "\n"; 
  }

  //FREE MEMORY
  delete new_rad_con;
}

//create diagrams of 1) Voronoi network and 2) accessible Voronoi network in an adjacent cell if specified with shift; xyz and vtk output files generated
void visVoro(char* name, double probeRad, int skel_a, int skel_b, int skel_c, VORONOI_NETWORK* vornet, ATOM_NETWORK* atmnet) {
  string filename_xyz = string(name).append("_voro.xyz");
  string filename2_xyz = string(name).append("_voro_accessible.xyz");
  string filename3_xyz = string(name).append("_voro_nonaccessible.xyz");
  vector<bool> accessInfo;
  vector<bool> nonaccessInfo;
  vector<CHANNEL> channels;
  CHANNEL::findChannels(vornet, probeRad, &accessInfo, &channels);
  int num_accessible = 0;
  for(unsigned int i = 0; i < accessInfo.size(); i++) if(accessInfo.at(i)==1) num_accessible++;

  int num_nonaccessible = 0;
  nonaccessInfo.resize(accessInfo.size(),false);
  for(unsigned int i = 0; i < accessInfo.size(); i++) 
    if(vornet->nodes.at(i).rad_stat_sphere>probeRad&&accessInfo.at(i)==0) 
       {
       num_nonaccessible++;
       nonaccessInfo[i]=true;
       };

  //write Voronoi nodes as .xyz
  FILE *output_xyz; output_xyz = fopen(filename_xyz.c_str(), "w");
  FILE *output2_xyz; output2_xyz = fopen(filename2_xyz.c_str(), "w");
  FILE *output3_xyz; output3_xyz = fopen(filename3_xyz.c_str(), "w");
  fprintf(output_xyz, "%d\nVoronoi diagram for %s with probe radius %.3f\n", (int)(accessInfo.size()), name, probeRad);
  fprintf(output2_xyz, "%d\nVoronoi accessible diagram for %s with probe radius %.3f\n", num_accessible, name, probeRad);
  fprintf(output3_xyz, "%d\nVoronoi non-accessible diagram for %s with probe radius %.3f\n", num_nonaccessible, name, probeRad);
  for(unsigned int i = 0; i < accessInfo.size(); i++) {
    fprintf(output_xyz, "Al %.3f %.3f %.3f %.3f\n", vornet->nodes.at(i).x, vornet->nodes.at(i).y, vornet->nodes.at(i).z, vornet->nodes.at(i).rad_stat_sphere); //"Al" for 'all' Voronoi nodes, shown in main unit cell
    if(accessInfo.at(i)==1) {
      Point abc = atmnet->xyz_to_abc(vornet->nodes.at(i).x, vornet->nodes.at(i).y, vornet->nodes.at(i).z);
      float new_a = abc.vals[0] + skel_a;
      float new_b = abc.vals[1] + skel_b;
      float new_c = abc.vals[2] + skel_c;
      Point new_xyz = atmnet->abc_to_xyz(new_a, new_b, new_c);
      fprintf(output2_xyz, "Ac %.3f %.3f %.3f %.3f\n", new_xyz.vals[0], new_xyz.vals[1], new_xyz.vals[2], vornet->nodes.at(i).rad_stat_sphere); //"Ac" for 'accessible' Voronoi nodes, shown in adjacent cell
    } else {
    if(vornet->nodes.at(i).rad_stat_sphere>probeRad) {
      Point abc = atmnet->xyz_to_abc(vornet->nodes.at(i).x, vornet->nodes.at(i).y, vornet->nodes.at(i).z);
      float new_a = abc.vals[0] + skel_a;
      float new_b = abc.vals[1] + skel_b;
      float new_c = abc.vals[2] + skel_c;
      Point new_xyz = atmnet->abc_to_xyz(new_a, new_b, new_c);
      fprintf(output3_xyz, "In %.3f %.3f %.3f %.3f\n", new_xyz.vals[0], new_xyz.vals[1], new_xyz.vals[2], vornet->nodes.at(i).rad_stat_sphere); //"In" for 'nonaccessible' Voronoi nodes, shown in adjacent cell

     }; // ends if for nodes larger than probeRad

    };
  }
  fclose(output_xyz);
  fclose(output2_xyz);
  fclose(output3_xyz);


  // VTK saving function
  // Note: to avoid issues with reindexing of nodes in .vtk file
  //       when printing (non)accesible network, all nodes are printed first

  //write Voronoi edges as .vtk
  string filename_vtk = string(name).append("_voro.vtk");
  string filename2_vtk = string(name).append("_voro_accessible.vtk");
  string filename3_vtk = string(name).append("_voro_nonaccessible.vtk");
  //start saving the addresses of (non)accessible nodes
  int *accessIndex;
  accessIndex = new int[accessInfo.size()];
  int *nonaccessIndex;
  nonaccessIndex = new int[nonaccessInfo.size()];
  int count = accessInfo.size(); //start counter at total number of nodes - because accessible nodes are added on after
  int countNonacc = nonaccessInfo.size(); //start counter at total number of nodes - because nonaccessible nodes are added on after
  for(unsigned int i = 0; i < accessInfo.size(); i++) {
    if(accessInfo.at(i)==1) {
      accessIndex[i] = count;
      count++;
    } else accessIndex[i] = -1;
    if(nonaccessInfo.at(i)==1) {
      nonaccessIndex[i] = countNonacc;
      countNonacc++;
    } else nonaccessIndex[i] = -1;
  }
  //refer to those addresses in edge storage
  vector <int> to_vector, from_vector;
  vector <int> to2_vector, from2_vector;
  vector <int> to3_vector, from3_vector;
  for(unsigned int i=0; i<vornet->edges.size(); i++) {
		if(accessInfo.at(vornet->edges.at(i).from)>0 && accessInfo.at(vornet->edges.at(i).to)>0 && vornet->edges.at(i).rad_moving_sphere>probeRad) {
      //this edge is accessible - save it for accessible part of the image
      if(vornet->edges.at(i).delta_uc_x==0 && vornet->edges.at(i).delta_uc_y==0 && vornet->edges.at(i).delta_uc_z==0) { //just edges within the unit cell
        from2_vector.push_back(accessIndex[vornet->edges.at(i).from]);
        to2_vector.push_back(accessIndex[vornet->edges.at(i).to]);
      } else if(vornet->edges.at(i).delta_uc_x==skel_a && vornet->edges.at(i).delta_uc_y==skel_b && vornet->edges.at(i).delta_uc_z==skel_c) { //just edges from overall section to accessible section
        from2_vector.push_back(vornet->edges.at(i).from);
        to2_vector.push_back(accessIndex[vornet->edges.at(i).to]);
      }
    } // ends if for accessible

      if(nonaccessInfo.at(vornet->edges.at(i).from)>0 && nonaccessInfo.at(vornet->edges.at(i).to)>0 && vornet->edges.at(i).rad_moving_sphere>probeRad) {
      //this edge is within nonaccessible pocket - save it for nonaccessible part of the image
      if(vornet->edges.at(i).delta_uc_x==0 && vornet->edges.at(i).delta_uc_y==0 && vornet->edges.at(i).delta_uc_z==0) { //just edges within the unit cell
        from3_vector.push_back(nonaccessIndex[vornet->edges.at(i).from]);
        to3_vector.push_back(nonaccessIndex[vornet->edges.at(i).to]);
      } else if(vornet->edges.at(i).delta_uc_x==skel_a && vornet->edges.at(i).delta_uc_y==skel_b && vornet->edges.at(i).delta_uc_z==skel_c) { //just edges from overall section to accessible section
        from3_vector.push_back(nonaccessIndex[vornet->edges.at(i).from]);
        to3_vector.push_back(nonaccessIndex[vornet->edges.at(i).to]);
      }
    } // ends if for nonaccessible

    if(vornet->edges.at(i).delta_uc_x==0 && vornet->edges.at(i).delta_uc_y==0 && vornet->edges.at(i).delta_uc_z==0) { //just edges within the unit cell
      //regardless of accessibility, want to save the within-cell edges for overall part of the image
      from_vector.push_back(vornet->edges.at(i).from);
      to_vector.push_back(vornet->edges.at(i).to);
    }
  }; // end loop over all edges


  // VTK file write
  FILE *output_vtk; output_vtk = fopen(filename_vtk.c_str(), "w");
  FILE *output2_vtk; output2_vtk = fopen(filename2_vtk.c_str(), "w");
  FILE *output3_vtk; output3_vtk = fopen(filename3_vtk.c_str(), "w");
  fprintf(output_vtk, "# vtk DataFile Version 2.0\nvtk data for file %s\nASCII\nDATASET POLYDATA\nPOINTS %d double\n", name, (int)(accessInfo.size()));
  fprintf(output2_vtk, "# vtk DataFile Version 2.0\nvtk data for file %s\nASCII\nDATASET POLYDATA\nPOINTS %d double\n", name, (int)(accessInfo.size())+num_accessible);
  fprintf(output3_vtk, "# vtk DataFile Version 2.0\nvtk data for file %s\nASCII\nDATASET POLYDATA\nPOINTS %d double\n", name, (int)(accessInfo.size())+num_nonaccessible);
  //write nodes to file
  for(unsigned int i = 0; i < accessInfo.size(); i++) {
    fprintf(output_vtk, "%.3f %.3f %.3f\n", vornet->nodes.at(i).x, vornet->nodes.at(i).y, vornet->nodes.at(i).z);
    fprintf(output2_vtk, "%.3f %.3f %.3f\n", vornet->nodes.at(i).x, vornet->nodes.at(i).y, vornet->nodes.at(i).z);
    fprintf(output3_vtk, "%.3f %.3f %.3f\n", vornet->nodes.at(i).x, vornet->nodes.at(i).y, vornet->nodes.at(i).z);
  }
  for(unsigned int i = 0; i < accessInfo.size(); i++) {
    if(accessInfo.at(i)==1) {
      Point abc = atmnet->xyz_to_abc(vornet->nodes.at(i).x, vornet->nodes.at(i).y, vornet->nodes.at(i).z);
      float new_a = abc.vals[0] + skel_a;
      float new_b = abc.vals[1] + skel_b;
      float new_c = abc.vals[2] + skel_c;
      Point new_xyz = atmnet->abc_to_xyz(new_a, new_b, new_c);
      fprintf(output2_vtk, "%.3f %.3f %.3f\n", new_xyz.vals[0], new_xyz.vals[1], new_xyz.vals[2]);
    };
    if(nonaccessInfo.at(i)==1) {
      Point abc = atmnet->xyz_to_abc(vornet->nodes.at(i).x, vornet->nodes.at(i).y, vornet->nodes.at(i).z);
      float new_a = abc.vals[0] + skel_a;
      float new_b = abc.vals[1] + skel_b;
      float new_c = abc.vals[2] + skel_c;
      Point new_xyz = atmnet->abc_to_xyz(new_a, new_b, new_c);
      fprintf(output3_vtk, "%.3f %.3f %.3f\n", new_xyz.vals[0], new_xyz.vals[1], new_xyz.vals[2]);
    };
  }
  //write edges to file
  fprintf(output_vtk, "LINES %d %d\n", (int)(to_vector.size()), ((int)(to_vector.size()))*3);
  for(unsigned int i = 0; i < to_vector.size(); i++) {
    fprintf(output_vtk, "2 %d %d\n", to_vector.at(i), from_vector.at(i));
  }
  fprintf(output2_vtk, "LINES %d %d\n", (int)(to2_vector.size()), ((int)(to2_vector.size()))*3);
  for(unsigned int i = 0; i < to2_vector.size(); i++) {
    fprintf(output2_vtk, "2 %d %d\n", to2_vector.at(i), from2_vector.at(i));
  }
  fprintf(output3_vtk, "LINES %d %d\n", (int)(to3_vector.size()), ((int)(to3_vector.size()))*3);
  for(unsigned int i = 0; i < to3_vector.size(); i++) {
    fprintf(output3_vtk, "2 %d %d\n", to3_vector.at(i), from3_vector.at(i));
  }
  fclose(output_vtk);
  fclose(output2_vtk);
  fclose(output3_vtk);
  //FREE MEMORY
  delete[] accessIndex;
}

//extract spherical substructures
void getLocalSubstructures(char* name, double probeRad, double local_substructure_radius, TRIPLET supercell_steps_local_substructure, string element, VORONOI_NETWORK* vornet, ATOM_NETWORK* atmnet, bool radial, bool usingElementOnly) {
  vector<bool> accessInfo;
  vector<CHANNEL> channels;
  CHANNEL::findChannels(vornet, probeRad, &accessInfo, &channels);
  ATOM_NETWORK simplified_atmnet;
  atmnet->copy(&simplified_atmnet);
  VORONOI_NETWORK simplified_vornet;
  vornet->copy(&simplified_vornet);
  vector<bool> isAccessible;
  int num_accessible = 0;
  //if not simplifying, avoid the calculation and just copy the data
  if(!usingElementOnly) {
    isAccessible = accessInfo;
  } else {
    //create the simplified atmnet, ready to calculate simplified vornet etc.
    simplified_vornet.edges.clear();
    simplified_vornet.nodes.clear();
    simplified_atmnet.atoms.clear();
    simplified_atmnet.numAtoms = 0;
    for(int i=0; i<atmnet->numAtoms; i++){
      if(atmnet->atoms.at(i).type.compare(element)==0) {
        ATOM newAtom;
        newAtom.type = atmnet->atoms.at(i).type;
        newAtom.a_coord = atmnet->atoms.at(i).a_coord;
        newAtom.b_coord = atmnet->atoms.at(i).b_coord;
        newAtom.c_coord = atmnet->atoms.at(i).c_coord;
        newAtom.x = atmnet->atoms.at(i).x;
        newAtom.y = atmnet->atoms.at(i).y;
        newAtom.z = atmnet->atoms.at(i).z;
        newAtom.radius = lookupRadius(newAtom.type, radial);
        simplified_atmnet.atoms.push_back(newAtom);
        simplified_atmnet.numAtoms++;
      }
    }
    cout << "Starting simplified Voronoi decomposition" << "\n";
    vector< VOR_CELL> simplified_cells;
    vector< BASIC_VCELL> simplified_bvcells;
    if(radial) {
      container_periodic_poly *rad_con = (container_periodic_poly *)performVoronoiDecomp(true, &simplified_atmnet, &simplified_vornet, simplified_cells, false, simplified_bvcells);
      delete rad_con;
    } else {
      container_periodic *no_rad_con = (container_periodic *)performVoronoiDecomp (false, &simplified_atmnet, &simplified_vornet, simplified_cells, false, simplified_bvcells); 
      delete no_rad_con;
    }
    cout << "Finished simplified Voronoi decomposition" << "\n";
    printf("*** Original voronoi network on %d atoms contains %d nodes and %d edges ***\n", atmnet->numAtoms, (int)(vornet->nodes.size()), (int)(vornet->edges.size()));
    printf("*** Simplified voronoi network on %d atoms contains %d nodes and %d edges ***\n", simplified_atmnet.numAtoms, (int)(simplified_vornet.nodes.size()), (int)(simplified_vornet.edges.size()));
    //now check which of these nodes of the simplified network lie in accessible space w.r.t. atmnet
    float chan_radius = probeRad;
    determineAccessibility(atmnet, chan_radius, probeRad, true, &isAccessible, &simplified_vornet); //where probeRad is the smaller, used to probe the channels accessible to chan_radius, the larger
  }
  for(int i=0; i<(int)(simplified_vornet.nodes.size()); i++) {
    if(simplified_vornet.nodes.at(i).rad_stat_sphere<=0) isAccessible.at(i)=false;
    if(isAccessible.at(i)==true) num_accessible++;
  }
  printf("*** %d of these %d nodes lie in accessible space ***\n", num_accessible, (int)(simplified_vornet.nodes.size()));
  //now try and remove redundant nodes
  int malloc_num_accessible = num_accessible;
  double **distance_between_nodes_array = new double*[malloc_num_accessible];
  for(int i=0; i<malloc_num_accessible; i++) distance_between_nodes_array[i] = new double[malloc_num_accessible];
  int *node_reference = new int[malloc_num_accessible];
  double *node_density = new double[malloc_num_accessible];
  int iter = 0;
  for(int i=0; i<(int)(simplified_vornet.nodes.size()); i++) {
    if(isAccessible.at(i)==true) {
      node_reference[iter] = i; //now we can iterate over this array for simplicity
      iter++;
    }
  }
  double minimum_node_distance = 1.0; //two nodes closer than this probably don't tell us anything different about the structure
  char some_distance_too_low_exists = 0;
  for(int i=0; i<malloc_num_accessible; i++) {
    node_density[i] = 0.0;
  }
  for(int i=0; i<malloc_num_accessible; i++) {
    distance_between_nodes_array[i][i] = 0.0;
    for(int j=i+1; j<malloc_num_accessible; j++) {
      int node_i = node_reference[i];
      int node_j = node_reference[j];
      Point abcCoord_i = atmnet->xyz_to_abc(simplified_vornet.nodes.at(node_i).x, simplified_vornet.nodes.at(node_i).y, simplified_vornet.nodes.at(node_i).z);
      Point abcCoord_j = atmnet->xyz_to_abc(simplified_vornet.nodes.at(node_j).x, simplified_vornet.nodes.at(node_j).y, simplified_vornet.nodes.at(node_j).z);
      double distance = atmnet->getDistCalc().minimum_periodic_distance(abcCoord_i.vals[0], abcCoord_i.vals[1], abcCoord_i.vals[2], abcCoord_j.vals[0], abcCoord_j.vals[1], abcCoord_j.vals[2]);
      if(distance<minimum_node_distance) some_distance_too_low_exists = 1;
      distance_between_nodes_array[i][j] = distance;
      distance_between_nodes_array[j][i] = distance;
      node_density[i] += exp(-1*distance*distance);
      node_density[j] += exp(-1*distance*distance);
    }
  }
  if(some_distance_too_low_exists>0) {
    //now we have the distance array, iterate over the nodes, removing the most dense, while there are some accessible ones closer than some distance
    while(some_distance_too_low_exists==1) {
      int accessibility_index_of_max_density_node = -1;
      double max_density = -1;
      for(int i=0; i<malloc_num_accessible; i++) {
        if(isAccessible.at(node_reference[i])==true) {
          if(node_density[i]>max_density) {
            max_density = node_density[i];
            accessibility_index_of_max_density_node = i;
          }
        }
      }
      if(max_density<0 || accessibility_index_of_max_density_node<0) {
        printf("ERROR: no node found in max density search\n");
        exit(EXIT_FAILURE);
      }
      //remove this most dense node
      isAccessible.at(node_reference[accessibility_index_of_max_density_node])=false;
      num_accessible--;
      //we need to recalculate densities and check if any remaining accessible node triggers the too-close test!
      some_distance_too_low_exists = 0;
      for(int i=0; i<malloc_num_accessible; i++) {
        node_density[i] = 0.0;
      }
      for(int i=0; i<malloc_num_accessible; i++) {
        int node_i = node_reference[i];
        if(isAccessible.at(node_i)==true) {
          for(int j=i+1; j<malloc_num_accessible; j++) {
            int node_j = node_reference[j];
            if(isAccessible.at(node_j)==true) {
              double distance = distance_between_nodes_array[i][j];
              if(distance<minimum_node_distance) some_distance_too_low_exists = 1;
              node_density[i] += exp(-1*distance*distance);
              node_density[j] += exp(-1*distance*distance);
            }
          }
        }
      }
    }
    printf("STATUS: %d nodes were removed, leaving %d, none of which are within %.3fA of any other\n", malloc_num_accessible-num_accessible, num_accessible, minimum_node_distance);
  }
  //free memory
  for(int i=0; i<malloc_num_accessible; i++) delete[] distance_between_nodes_array[i];
  delete[] distance_between_nodes_array;
  delete[] node_reference;
  delete[] node_density;

  //construct the supercell_atmnet for extracting spheres
  ATOM_NETWORK supercell_atmnet;
  atmnet->copy(&supercell_atmnet);
  for(int a_step=-1*supercell_steps_local_substructure[0]; a_step<=supercell_steps_local_substructure[0]; a_step++) {
    for(int b_step=-1*supercell_steps_local_substructure[1]; b_step<=supercell_steps_local_substructure[1]; b_step++) {
      for(int c_step=-1*supercell_steps_local_substructure[2]; c_step<=supercell_steps_local_substructure[2]; c_step++) {
        if(a_step!=0 || b_step!=0 || c_step!=0) { //don't make copies within the unit cell!
          for(int i=0; i<atmnet->numAtoms; i++){
            ATOM newAtom;
            newAtom.type = atmnet->atoms.at(i).type;
            newAtom.a_coord = atmnet->atoms.at(i).a_coord+a_step;
            newAtom.b_coord = atmnet->atoms.at(i).b_coord+b_step;
            newAtom.c_coord = atmnet->atoms.at(i).c_coord+c_step;
            Point xyz = atmnet->abc_to_xyz(newAtom.a_coord, newAtom.b_coord, newAtom.c_coord);
            newAtom.x = xyz.vals[0];
            newAtom.y = xyz.vals[1];
            newAtom.z = xyz.vals[2];
            newAtom.radius = lookupRadius(newAtom.type, radial);
            supercell_atmnet.atoms.push_back(newAtom);
            supercell_atmnet.numAtoms++;
          }
        }
      }
    }
  }
  printf("*** Unit cell contains %d atoms; supercell contains %d atoms ***\n", atmnet->numAtoms, supercell_atmnet.numAtoms);
  int node_index = 0;
  for(int i=0; i<(int)(simplified_vornet.nodes.size()); i++) {
    vector<ATOM> vector_of_local_atoms; //all those visible atoms within some distance, i.e. those which comprise the local surface
    if(isAccessible.at(i)==true) {
      for(int j=0; j<supercell_atmnet.numAtoms; j++) {
        float distance = calcEuclideanDistance(simplified_vornet.nodes.at(i).x, simplified_vornet.nodes.at(i).y, simplified_vornet.nodes.at(i).z, supercell_atmnet.atoms.at(j).x, supercell_atmnet.atoms.at(j).y, supercell_atmnet.atoms.at(j).z);
        if(distance<=local_substructure_radius) vector_of_local_atoms.push_back(supercell_atmnet.atoms.at(j));
      }
      //now finished putting local atoms in the vector, write them to file for now
      string number = intAsString(node_index);
      string filename_xyz = string(name).append("_local_substructure_");
      filename_xyz.append(number);
      filename_xyz.append(".xyz");
      FILE *local_xyz; local_xyz = fopen(filename_xyz.c_str(), "w");
      fprintf(local_xyz, "%d\nxyz header\n", (int)vector_of_local_atoms.size()+1);
      for(int j=0; j<(int)vector_of_local_atoms.size(); j++) {
        fprintf(local_xyz, "%s %.3f %.3f %.3f\n", vector_of_local_atoms.at(j).type.c_str(), vector_of_local_atoms.at(j).x, vector_of_local_atoms.at(j).y, vector_of_local_atoms.at(j).z);
      }
      fprintf(local_xyz, "X %.3f %.3f %.3f\n", simplified_vornet.nodes.at(i).x, simplified_vornet.nodes.at(i).y, simplified_vornet.nodes.at(i).z);
      fclose(local_xyz);
      node_index++;
    }
  }
}


double calcASA(ATOM_NETWORK *hiaccatmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double r_probe_chan, double r_probe, double rho_crystal, int numSamples, bool excludePockets, ostream &output, char *filename, bool visualize, bool VisITflag, bool LiverpoolFlag, bool ExtendedOutputFlag){

  ATOM_NETWORK *atmnet; // pointed to the analyzed (original atomsnet)

  if(highAccuracy) atmnet=orgatmnet; else atmnet = hiaccatmnet; // this ensures that the loops below loop over original atoms
                                                                // and not "high-acuracy" clusters, which are used only for
                                                                // voronoi decomposition

  // Create an object that handles analysis of accessibility of sampled points
  AccessibilityClass accessAnalysis;
  if(highAccuracy) accessAnalysis.setupAndFindChannels(hiaccatmnet, orgatmnet, highAccuracy, r_probe_chan, r_probe);
    else accessAnalysis.setupAndFindChannels(hiaccatmnet, hiaccatmnet, highAccuracy, r_probe_chan, r_probe);

  accessAnalysis.removeOverlappedNodes();

  srand(randSeed);

  vector<Point> axsPoints   = vector<Point> (); // List of accessible points to be visualized in ZeoVis
  vector<int> axsPointsChannelIDs;  // corresponding channel IDs
  vector<Point> inaxsPoints = vector<Point> (); // List of inaccessible points to be visualized in ZeoVis
  vector<int> inaxsPointsPocketIDs; // corresponding pocket IDs
  double totalSA = 0;
  double totalSA_inaxs = 0;

  vector<double> totalSA_inChannel(accessAnalysis.n_channels,0);
  vector<double> totalSA_inPocket(accessAnalysis.n_pockets,0);

  //vectors that store SA contributions per atom (a pair of accessible and inaccessible)
  vector< pair <double,double> > histogramSAperAtom;

  // Sample around each of the network's atoms
  for(int i = 0; i < atmnet->numAtoms; i++){
    int count = 0;
    int count_inaxs = 0;
    // MC statistics can be collected w.r.t to identified channels and inaccessible pockets
    vector<int> count_inChannel(accessAnalysis.n_channels,0);
    vector<int> count_inPocket(accessAnalysis.n_pockets,0);

    for(int j = 0; j < numSamples; j++){
      bool overlaps = false;
      bool inside = false;

      // Randomly sample point on a sphere of radius 1
      double theta = (rand()*1.0/RAND_MAX)*2*PI;
      double cosphi = 1.0 - (rand()*1.0/RAND_MAX)*2.0;
      double phi = acos(cosphi);

      // Convert spherical coordinates to xyz coordinates
      double xpoint = sin(phi)*cos(theta); //Don't need abs(sin(phi)) becase phi is from 0 to PI
      double ypoint = sin(phi)*sin(theta);
      double zpoint = cosphi;	
      
      // Scale the point to lie on a radius of r_probe + r_atom#i
      xpoint *= (atmnet->atoms[i].radius + r_probe);
      ypoint *= (atmnet->atoms[i].radius + r_probe);
      zpoint *= (atmnet->atoms[i].radius + r_probe);

      // Convert (x,y,z) coordinates to ones relative to the unit cell vectors
      Point abc_coords = atmnet->xyz_to_abc(xpoint, ypoint, zpoint);

      // Transform coords relative to the atom center into ones relative to the unit cell vectors
      double newAPoint = abc_coords[0] + atmnet->atoms[i].a_coord;
      double newBPoint = abc_coords[1] + atmnet->atoms[i].b_coord;
      double newCPoint = abc_coords[2] + atmnet->atoms[i].c_coord;

      Point samplingPoint = atmnet->abc_to_xyz(newAPoint, newBPoint, newCPoint);


      // Calling accessibility object to determine accessibility of the point (this replaced a big chunk of code by Thomas)
      pair<bool,bool> answer = (accessAnalysis.isSPointInsideAtomAndNotAccessible(samplingPoint, i));
      inside = answer.first; overlaps = answer.second;
      if(accessAnalysis.needToResample() == true) j--; // the sampled point could not be analyzed in isSPointInsideAtomAndNotAccessible() function, resampling needed

      if(inside == false && excludePockets == false) overlaps = false; // if ignore inacceible pockets, treat the point as accessible (unless inside atom)


      // Store sampled points that did not overlap within atom for later visualization
      if(accessAnalysis.needToResample() == false && !inside) {
	Point abcCoords = Point(newAPoint, newBPoint, newCPoint);
	Point coords = atmnet->abc_to_xyz(atmnet->shiftABCInUC(abcCoords));
	if(overlaps)
          {
          count_inaxs++;
	  inaxsPoints.push_back(coords);
          pair <int,int> CoP = accessAnalysis.lastChannelOrPocket();
          if(CoP.first!=-1)
            {
            cout << "Error: CoP.first!=-1 in pocket, consult source code provider\n";
            } else {
            count_inPocket[CoP.second]++;
            };
          inaxsPointsPocketIDs.push_back(CoP.second);
          }
	else
          {
	  axsPoints.push_back(coords);
          count++;
          pair <int,int> CoP = accessAnalysis.lastChannelOrPocket();
          if(CoP.second!=-1)
            {
            cout << "Error: CoP.second!=-1 in channel, consult source code provider\n";
            } else {
            count_inChannel[CoP.first]++;
            };
          axsPointsChannelIDs.push_back(CoP.first);
          };
      }

//      if(accessAnalysis.needToResample() == false && !overlaps)
//        count++;
    }
   
    // SA per atom (in accessible space and in inaccessible pockets)
    double SAperAtom = (1.0*count)/numSamples * 4.0 * PI * pow((atmnet->atoms.at(i).radius + r_probe),2);
    double SAperAtom_inaxs = (1.0*count_inaxs)/numSamples * 4.0 * PI * pow((atmnet->atoms.at(i).radius + r_probe),2);

    //Increment total surface area by fraction of atom that is accessible
    totalSA += SAperAtom;
    //and inaccessible
    totalSA_inaxs += SAperAtom_inaxs;

    // SA per channel and pocket
    for(unsigned int m = 0; m < totalSA_inChannel.size(); m++)
      {
      totalSA_inChannel[m] += (1.0*count_inChannel[m])/numSamples * 4.0 * PI * pow((atmnet->atoms.at(i).radius + r_probe),2);
      };
    for(unsigned int m = 0; m < totalSA_inPocket.size(); m++)
      {
      totalSA_inPocket[m] += (1.0*count_inPocket[m])/numSamples * 4.0 * PI * pow((atmnet->atoms.at(i).radius + r_probe),2);
      };

    histogramSAperAtom.push_back(pair <double,double> (SAperAtom,SAperAtom_inaxs));

  }

  // Write the necessary commands to the output stream
  // necessary to visualize the sampling
  if(visualize){
    if(VisITflag == false)
      {
      reportPoints(output, axsPoints, inaxsPoints);
      }
    else{
      if(LiverpoolFlag == false)
        {
        //report points in std. VisIt format
        reportPointsVisIT(output, axsPoints, inaxsPoints);
        }
      else{
      //report points in Liverpool format
      vector<Point> LiverpoolAxsPoints = vector<Point> ();
      vector<Point> LiverpoolInaxsPoints = vector<Point> ();
      for(unsigned i = 0; i < axsPoints.size(); i++)
         {
         LiverpoolAxsPoints.push_back(atmnet->xyz_to_abc(axsPoints[i]));
         };
      for(unsigned j = 0; j < inaxsPoints.size(); j++)
         {
         LiverpoolInaxsPoints.push_back(atmnet->xyz_to_abc(inaxsPoints[j]));
         };
//      reportPointsVisIT(output, LiverpoolAxsPoints, LiverpoolInaxsPoints);
      reportPointsVisIT(output, LiverpoolAxsPoints, axsPointsChannelIDs, LiverpoolInaxsPoints, inaxsPointsPocketIDs);
      };
      };
    //reportResampledPoints(output, resampledInfo);
  }


  double ucVolume    = calcDeterminant(atmnet->ucVectors);   
  double saPerVolume = (totalSA/ucVolume)*pow(10.0,4.0); // Units of m^2/cm^3
  double saPerMass   = saPerVolume/rho_crystal;          // Units of m^2/g

  double saPerVolume_inaxs = (totalSA_inaxs/ucVolume)*pow(10.0,4.0); // Units of m^2/cm^3
  double saPerMass_inaxs   = saPerVolume_inaxs/rho_crystal;          // Units of m^2/g

  if(!visualize){
//    output << newAtomNet.name << " ";
    output << "@ " << filename << " ";
    output << "Unitcell_volume: " << ucVolume << "   Density: " << rho_crystal << "   ";
    output << "ASA_A^2: " << totalSA << " "
	   << "ASA_m^2/cm^3: " << saPerVolume << " "
	   << "ASA_m^2/g: " << saPerMass << " "
           << "NASA_A^2: " << totalSA_inaxs << " "
           << "NASA_m^2/cm^3: " << saPerVolume_inaxs << " "
           << "NASA_m^2/g: " << saPerMass_inaxs << "\n";
  }

  // Surface composition histogram
  if(!visualize&&ExtendedOutputFlag){
    double metalFrac = 0, metalFrac_inaxs = 0;

    for(int i = 0; i < atmnet->numAtoms; i++)
       {
       if(isMetal(atmnet->atoms[i].type) == true) 
         {
         metalFrac+= histogramSAperAtom.at(i).first;
         metalFrac_inaxs += histogramSAperAtom.at(i).second;
         };
       };

    metalFrac = metalFrac / totalSA;
    metalFrac_inaxs = metalFrac_inaxs / totalSA_inaxs;

    output << "Metal fraction in ASA and NASA: " << metalFrac << "   " << metalFrac_inaxs << "\n";
  };

  // extended output
  if(!visualize){
    output << "Number_of_channels: " << totalSA_inChannel.size()  << " Channel_surface_area_A^2: ";
    for(unsigned int i = 0; i < totalSA_inChannel.size(); i++)
      {
      output << totalSA_inChannel[i] << "  ";
      };
    output << "\nNumber_of_pockets: " << totalSA_inPocket.size()  << " Pocket_surface_area_A^2: ";
    for(unsigned int i = 0; i < totalSA_inPocket.size(); i++)
      {
      output << totalSA_inPocket[i] << "  ";
      };
    output << "\n";
    };

    
    // Warn user if points were resampled
  int resampleCount =  accessAnalysis.getResampleCount();
  if(resampleCount != 0){
    cerr << "\n" << "\n"
  	 << "Warning: Resampled " << resampleCount << " points out of " << atmnet->numAtoms*numSamples 
	 << "\n" << "\n"; 
  }
  accessAnalysis.deconstruct();
  return totalSA;
}


/* backup of ASA function


double calcASA(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double r_probe_chan, double r_probe, double rho_crystal, int numSamples, bool excludePockets, ostream &output, char *filename, bool visualize, bool VisITflag, bool LiverpoolFlag){

  // Create a temporary copy of the atomic network in which each atom's radius has been increased by the probe radius
  ATOM_NETWORK newAtomNet;
  atmnet->copy(&newAtomNet);
  for(int i = 0; i < newAtomNet.numAtoms; i++){ newAtomNet.atoms[i].radius += r_probe; }  

  // Calculate and store the Voronoi network for this new atomic network
  VORONOI_NETWORK vornet;
  vector<BASIC_VCELL> vorcells;
  vector<VOR_CELL> advCells;
  container_periodic_poly *new_rad_con = (container_periodic_poly *)performVoronoiDecomp(true, &newAtomNet, &vornet, advCells, false, vorcells);

  vector<CHANNEL> channels = vector<CHANNEL>();
  vector<bool> accessInfo = vector<bool> ();
  CHANNEL::findChannels(&vornet, max(0.0, r_probe_chan - r_probe), &accessInfo, &channels);
  srand(randSeed);

  // Remove all nodes from Voronoi cells that lie within the sampling sphere
  for(unsigned int i = 0; i < vorcells.size(); i++){
    vorcells[i].removeOverlappedNodes(i, atmnet, r_probe);
  }

  vector<Point> axsPoints   = vector<Point> (); // List of accessible points to be visualized in ZeoVis
  vector<Point> inaxsPoints = vector<Point> (); // List of inaccessible points to be visualized in ZeoVis
  vector< pair<int, Point> > resampledInfo = vector< pair<int, Point> > ();  // List of resampled points and the id of the Voronoi cell to which they belong
  double totalSA = 0;
  double totalSA_inaxs = 0;
  int resampleCount = 0;

  //vectors that store SA contributions per atom (a pair of accessible and inaccessible)
  vector< pair <double,double> > histogramSAperAtom;

  // Sample around each of the network's atoms
  for(int i = 0; i < atmnet->numAtoms; i++){
    int count = 0;
    int count_inaxs = 0;
    for(int j = 0; j < numSamples; j++){
      bool overlaps = false;

      // Randomly sample point on a sphere of radius 1
      double theta = (rand()*1.0/RAND_MAX)*2*PI;
      double cosphi = 1.0 - (rand()*1.0/RAND_MAX)*2.0;
      double phi = acos(cosphi);

      // Convert spherical coordinates to xyz coordinates
      double xpoint = sin(phi)*cos(theta); //Don't need abs(sin(phi)) becase phi is from 0 to PI
      double ypoint = sin(phi)*sin(theta);
      double zpoint = cosphi;	
      
      // Scale the point to lie on a radius of r_probe + r_atom#i
      xpoint *= (atmnet->atoms[i].radius + r_probe);
      ypoint *= (atmnet->atoms[i].radius + r_probe);
      zpoint *= (atmnet->atoms[i].radius + r_probe);

      // Convert (x,y,z) coordinates to ones relative to the unit cell vectors
      Point abc_coords = atmnet->xyz_to_abc(xpoint, ypoint, zpoint);

      // Transform coords relative to the atom center into ones relative to the unit cell vectors
      double newAPoint = abc_coords[0] + atmnet->atoms[i].a_coord;
      double newBPoint = abc_coords[1] + atmnet->atoms[i].b_coord;
      double newCPoint = abc_coords[2] + atmnet->atoms[i].c_coord;

      Point samplingPoint = atmnet->abc_to_xyz(newAPoint, newBPoint, newCPoint);

      double atomX, atomY, atomZ; // Coordinates of atom closest to point as calculated in the next line

      int minID; 
      bool foundCell = new_rad_con->find_voronoi_cell(samplingPoint[0], samplingPoint[1], samplingPoint[2], atomX, atomY, atomZ, minID);
      if(!foundCell){
	cerr << "Error: Unable to find Voronoi cell for sampled point in ASA calculation." << "\n"
	     << "Occurred for structure " << newAtomNet.name << "\n"
	     << "Exiting..." << "\n";
	exit(1);
      }

      // If point in Voronoi cell of different atom, probe-atom overlap occurs because of d^2-r^2 criterion.
      if(minID != i)
	overlaps = true;

      // Check for overlap with periodic image of original sampling sphere
      ATOM curAtom = atmnet->atoms[i];

      // Check for overlap with periodic image of original sampling sphere
      double dx = samplingPoint[0] - atomX;
      double dy = samplingPoint[1] - atomY;
      double dz = samplingPoint[2] - atomZ;
      double minDist = sqrt(dx*dx + dy*dy + dz*dz);
      if(minDist < curAtom.radius + r_probe - 0.00000001){
	overlaps = true;
      }

bool inside = overlaps; // Only want to visualize points that are not inside sampling spheres

      // If no probe-atom overlap and pockets are being excluded, check Voronoi nodes of cell to determine accessibility of point
      if(!overlaps && excludePockets){
	Point circCenter   = Point(curAtom.x, curAtom.y, curAtom.z);
	Point sampleRay    = Point(xpoint, ypoint, zpoint);
                                                                                
	// Select the relevant Voronoi cell
	bool foundNode = false;
	BASIC_VCELL vcell = vorcells[i];

	if(vcell.getNumNodes() == 0){
	  cerr << "Error occurred during ASA calculations: Voronoi cell of sampled point does not have any nodes." << "\n"
	       << "Voronoi cell is #" << minID << " in structure " << newAtomNet.name << "\n"
	       << "Please contact the source code provider with this input." << "\n"
	       << "Exiting..." << "\n";
	  exit(1);
	}
     
	// Scan the nodes in the Voronoi cell to find if a line can be drawn from the node to the sampling point
	for(int k = 0; k < vcell.getNumNodes(); k++){
	  Point nodePoint = vcell.getNodeCoord(k);
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
	    
	    // Access status of node determines sampled point's accessibility
	    overlaps = !accessInfo.at(vcell.getNodeID(k));
	    break;
	  }
	}

	// Sampling failed due to lying on Voronoi cell face and numerical inaccurarcy. 
	// Record failure, resample and notify user after calculation is finished
	if(!foundNode){
	  resampleCount++;
	  resampledInfo.push_back(pair<int, Point> (minID, samplingPoint));
	  j -= 1;
	}
      }

      // Store sampled points that did not overlap within atom for later visualization
      if(!inside) {
	Point abcCoords = Point(newAPoint, newBPoint, newCPoint);
	Point coords = atmnet->abc_to_xyz(atmnet->shiftABCInUC(abcCoords));
	if(overlaps)
          {
          count_inaxs++;
	  inaxsPoints.push_back(coords);
          }
	else
	  axsPoints.push_back(coords);
      }

      if(!overlaps)
        count++;
    }
   
    // SA per atom (in accessible space and in inaccessible pockets)
    double SAperAtom = (1.0*count)/numSamples * 4.0 * PI * pow((atmnet->atoms.at(i).radius + r_probe),2);
    double SAperAtom_inaxs = (1.0*count_inaxs)/numSamples * 4.0 * PI * pow((atmnet->atoms.at(i).radius + r_probe),2);

    //Increment total surface area by fraction of atom that is accessible
    totalSA += SAperAtom;
    //and inaccessible
    totalSA_inaxs += SAperAtom_inaxs;

    histogramSAperAtom.push_back(pair <double,double> (SAperAtom,SAperAtom_inaxs));

  }

  // Write the necessary commands to the output stream
  // necessary to visualize the sampling
  if(visualize){
    if(VisITflag == false)
      {
      reportPoints(output, axsPoints, inaxsPoints);
      }
    else{
      if(LiverpoolFlag == false)
        {
        //report points in std. VisIt format
        reportPointsVisIT(output, axsPoints, inaxsPoints);
        }
      else{
      //report points in Liverpool format
      vector<Point> LiverpoolAxsPoints = vector<Point> ();
      vector<Point> LiverpoolInaxsPoints = vector<Point> ();
      for(unsigned i = 0; i < axsPoints.size(); i++)
         {
         LiverpoolAxsPoints.push_back(atmnet->xyz_to_abc(axsPoints[i]));
         };
      for(unsigned j = 0; j < inaxsPoints.size(); j++)
         {
         LiverpoolInaxsPoints.push_back(atmnet->xyz_to_abc(inaxsPoints[j]));
         };
      reportPointsVisIT(output, LiverpoolAxsPoints, LiverpoolInaxsPoints);
      };
      };
    //reportResampledPoints(output, resampledInfo);
  }


  double ucVolume    = calcDeterminant(atmnet->ucVectors);   
  double saPerVolume = (totalSA/ucVolume)*pow(10.0,4.0); // Units of m^2/cm^3
  double saPerMass   = saPerVolume/rho_crystal;          // Units of m^2/g

  double saPerVolume_inaxs = (totalSA_inaxs/ucVolume)*pow(10.0,4.0); // Units of m^2/cm^3
  double saPerMass_inaxs   = saPerVolume_inaxs/rho_crystal;          // Units of m^2/g

  if(!visualize){
//    output << newAtomNet.name << " ";
    output << "@ " << filename << " ";
    output << "ASA_A^2: " << totalSA << " "
	   << "ASA_m^2/cm^3: " << saPerVolume << " "
	   << "ASA_m^2/g: " << saPerMass << " "
           << "NASA_A^2: " << totalSA_inaxs << " "
           << "NASA_m^2/cm^3: " << saPerVolume_inaxs << " "
           << "NASA_m^2/g: " << saPerMass_inaxs << "\n";
  }

  // Surface composition histogram
  if(!visualize){
    double metalFrac = 0, metalFrac_inaxs = 0;

    for(int i = 0; i < atmnet->numAtoms; i++)
       {
       if(isMetal(atmnet->atoms[i].type) == true) 
         {
         metalFrac+= histogramSAperAtom.at(i).first;
         metalFrac_inaxs += histogramSAperAtom.at(i).second;
         };
       };

    metalFrac = metalFrac / totalSA;
    metalFrac_inaxs = metalFrac_inaxs / totalSA_inaxs;

    output << "Metal fraction in ASA and NASA: " << metalFrac << "   " << metalFrac_inaxs << "/n";
  };
    
    // Warn user if points were resampled
  if(resampleCount != 0){
    cerr << "\n" << "\n"
  	 << "Warning: Resampled " << resampleCount << " points out of " << atmnet->numAtoms*numSamples 
	 << "\n" << "\n"; 
  }
  delete new_rad_con;
  return totalSA;
}


// ends of asa backup
*/


double calcASA(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double r_probe, double rho_crystal, int numSamples, bool excludePockets, ostream &output, char *filename, bool visualize, bool VisITflag, bool LiverpoolFlag, bool ExtendedOutputFlag){
  return calcASA(atmnet, orgatmnet, highAccuracy, r_probe, r_probe, rho_crystal, numSamples, excludePockets, output, filename, visualize, VisITflag, LiverpoolFlag, ExtendedOutputFlag);
}

/// Cython wrapper to calcASA where visualization flags are ignored
string calcASA(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double r_probe_chan, double r_probe, int numSamples, bool excludePockets, bool ExtendedOutputFlag)
{
    stringstream output;
    string filename = "No filename";
    double rho_crystal = calcDensity(atmnet);
    double sa = calcASA(atmnet, orgatmnet, highAccuracy, r_probe_chan, r_probe, rho_crystal, numSamples, excludePockets, output, (char *)filename.data(), false, false, false, ExtendedOutputFlag);
    return output.str();
}


/** Returns the volume accessible to a particle of the provided radius. Accessible volume is defined as any region of space
 *  in which the center of the particle of the provided radius can reach. Excludes inaccessible pockets if requested .*/
double calcAV(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double r_probe, int numSamples, bool excludePockets, ostream &output, char *filename, bool visualize, bool VisITflag, bool LiverpoolFlag, bool blockingMode, double low_dist_cutoff, double high_dist_cutoff){
  return calcAV(atmnet, orgatmnet, highAccuracy, r_probe, r_probe, numSamples, excludePockets, output, filename, visualize, VisITflag, LiverpoolFlag, blockingMode, low_dist_cutoff, high_dist_cutoff);
}

/// Cython wrapper to calcAV where visualization flags are ignored
string calcAV(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double r_probe_chan, double r_probe, int numSamples, bool excludePockets,  double low_dist_cutoff, double high_dist_cutoff)
{
    stringstream output;
    string filename = "No filename";
    double vol = calcAV(atmnet, orgatmnet, highAccuracy, r_probe_chan, r_probe, numSamples, excludePockets, output, (char *)filename.data(), false, false, false, false, low_dist_cutoff, high_dist_cutoff);
    return output.str();
}



/* accessiblePoint function - to be removed (backup only) */

/** Returns if a given point is accessible. Accessible points are defined as points that lie inside of the void network that is accessible to a probe of radius r_probe.*/
/*  This function was added by Christopher based on CalcAV function provided by Thomas; to be removed after adding accessibility class */
/*
bool accessiblePoint(Point sample_pnt,container_periodic_poly *new_rad_con, double r_probe, ATOM_NETWORK *atmnet, vector<BASIC_VCELL> vorcells,vector<bool> accessInfo){
  bool overlaps = false;
  bool excludePockets = true;
  double newAtomX, newAtomY, newAtomZ;
  int minAtomID;
  bool foundCell = new_rad_con->find_voronoi_cell(sample_pnt[0],sample_pnt[1],sample_pnt[2], newAtomX, newAtomY, newAtomZ, minAtomID);
  if(!foundCell){
    cerr << "Error: Unable to find Voronoi cell for sampled point in AV calculation." << "\n"
	 << "Occurred for structure " << atmnet->name << "\n"
	 << "Exiting..." << "\n";
    exit(1);
  }
  
  ATOM curAtom = atmnet->atoms[minAtomID];
  
  // Adjust sampling point so that it lies within the Voronoi cell of interest constructed previously
  sample_pnt = (sample_pnt.add(Point(curAtom.x, curAtom.y, curAtom.z).subtract(Point(newAtomX, newAtomY, newAtomZ))));
  
  double minDist = calcEuclideanDistance(sample_pnt[0], sample_pnt[1], sample_pnt[2], curAtom.x, curAtom.y, curAtom.z);
  
  if(minDist < r_probe + curAtom.radius - threshold)
    overlaps = true;
  
  bool inside = overlaps;
  
  // If necessary, check Voronoi nodes of cell to determine accessibility of point
  if(!overlaps && excludePockets){
    BASIC_VCELL vcell = vorcells[minAtomID];
    Point circCenter = Point(curAtom.x, curAtom.y, curAtom.z);
    double samplingRadius = minDist;
    Point samplePoint = Point(sample_pnt[0]-curAtom.x, sample_pnt[1]-curAtom.y, sample_pnt[2]-curAtom.z);
    
    // Scan the nodes in the Voronoi cell to find if line can be drawn from the node to the sampling point
    bool foundNode = false;
    if(vcell.getNumNodes() == 0){
      cerr << "Error: Voronoi cell of sampled point does not have any nodes" << "\n"
	   << "Point: " << sample_pnt << "\n"
	   << "Voronoi cell is #" << minAtomID << " in structure " << atmnet->name << "\n"
	   << "Please contact the source code provider." << "\n"
	   << "Exiting..." << "\n";
      exit(1);
    }
    for(int k = 0; k < vcell.getNumNodes(); k++){
      Point nodePoint = vcell.getNodeCoord(k);  
      bool nodeInsideSphere = (calcEuclideanDistance(nodePoint, circCenter) < samplingRadius);
      if(!nodeInsideSphere){
	if(samplePoint*(sample_pnt-nodePoint) > 0) {
	  // Angle is less than 90 degrees and so the line segment intersects twice,
	  // making the path not viable
	}
	else {
	  // Angle is at least 90 degrees and so the line segment interesects only once, 
	  // thereby representing a viable path
	  foundNode = true;
	  overlaps = !accessInfo.at(vcell.getNodeID(k));
	  break;
	}
      }
    } 
    if (foundNode == true){
      return !overlaps;
    }
    else {
      cerr << "Error: No node found for associted point in Vorcell" << sample_pnt << endl;
    }  
  } 
  return false;
}
*/




/** Blocking function 
    Given a set of accessible and inaccessible points, generate blocking spheres and save to output stream
    Input arguments have vectors: coordinates of points (separate accesible and inaccessible) and vectors with pore IDs (to idnetify to which channel/pocket a point belongs)
    currently points are in factional coordinates
**/

void blockPockets(ATOM_NETWORK *atmnet, ostream &output, vector<Point> axsPoints, vector<int> axsPChIDs, vector<Point> inaxsPoints, vector<int> inaxsPPIDs, double probeRad){
/*
  for(unsigned int i = 0; i < axsPoints.size(); i++){
    Point coords = axsPoints.at(i);
    output << coords[0] << " " << coords[1] << " " << coords[2] << " 1 a " << axsPChIDs[i] << "\n";
   }
*/
/*
  for(unsigned int i = 0; i < inaxsPoints.size(); i++){
    Point coords = inaxsPoints.at(i);
    // print coordinates and pocket ID
    output <<  coords[0] << " " << coords[1] << " " << coords[2] << " 0 block " << inaxsPPIDs[i] << "\n";
  }
*/

  //DESCRIPTION: iterate over each pore ID that corresponds to a pocket, and block all points in that pocket with spheres without blocking any non-pocket points - the spheres are created iteratively, and many spheres may be needed per pocket
  bool debug = false; //prints extra info to terminal, including .xyz format Cartesian pocket info for visualisation
  vector<SPHERE> spheres_vector; //stores blocking spheres we have created below;
  double sphere_radius_overshoot = 0.1; //optional - the amount in A by which a sphere can extend beyond the furthest point to be blocked (so that the sphere does not perfectly intersect the furthest MC point, potentially leaving volume nearby that should have been blocked)

  //1) identify how many pores we have, and which ones are pockets
  int max_pore_ID = 0;
  int num_axsPoints = axsPoints.size();
  int num_inaxsPoints = inaxsPoints.size();
  for(int i=0; i<num_axsPoints; i++) if(axsPChIDs.at(i)>max_pore_ID) max_pore_ID = axsPChIDs.at(i);
  for(int i=0; i<num_inaxsPoints; i++) if(inaxsPPIDs.at(i)>max_pore_ID) max_pore_ID = inaxsPPIDs.at(i);
  int num_pores = max_pore_ID+1;
  vector<bool> is_pocket;
  for(int i=0; i<num_pores; i++) is_pocket.push_back(false);
  for(int i=0; i<num_inaxsPoints; i++) if(!is_pocket.at(inaxsPPIDs.at(i))) is_pocket.at(inaxsPPIDs.at(i)) = true;
  if(debug) {
    printf("DEBUG: there are %d pores and they are assigned to channels and pockets as follows:\n", num_pores);
    for(int i=0; i<num_pores; i++) {
      printf("ID %d: ", i);
      if(is_pocket.at(i)) printf("pocket\n"); else printf("channel\n");
    }
  }

  //2) for each pocket, determine which points are members of this pocket and not already blocked by an earlier sphere (hence the resulting points may not all be connected), and store them all in selected_points
  for(int i=0; i<num_pores; i++) {
    vector<Point> selected_points; //stores to-be-blocked points in this pore
    if(is_pocket.at(i)) { //this pore is a pocket - we are going to block it with spheres
      for(int j=0; j<num_inaxsPoints; j++) {
        if(inaxsPPIDs.at(j)==i) { //inaxs point j is in pore i, so needs to be blocked - if it is not already blocked, push to selected_points vector
          Point p = inaxsPoints.at(j);
          bool already_blocked = false;
          int num_spheres = spheres_vector.size();
					for(int s=0; s<num_spheres && !already_blocked; s++) {
            SPHERE this_sphere = spheres_vector.at(s);
            double dist = atmnet->calcDistanceABC(p[0], p[1], p[2], this_sphere.x, this_sphere.y, this_sphere.z);
						if(dist<this_sphere.r) already_blocked = true;
					}
          if(!already_blocked) selected_points.push_back(p);
        }
      }

      //3) begin a loop - while there are points in selected_points, we need to do more blocking
      int num_selected_points = selected_points.size();
      while(num_selected_points>0) {
        if(debug) printf("DEBUG: there are %d points left to be blocked in pore with ID %d\n", num_selected_points, i);
        int most_dense_index = get_most_dense_index(atmnet, &selected_points);
        Point most_dense = selected_points.at(most_dense_index); //find most dense point in the remaining point cloud
        double closest_channel = -1; //distance to nearest MC point that is in a channel
        for(int j=0; j<num_axsPoints; j++) { //loop over points to be preserved, and find closest
          Point p = axsPoints.at(j);
          double dist = atmnet->calcDistanceABC(p[0], p[1], p[2], most_dense[0], most_dense[1], most_dense[2]);
          if(dist<closest_channel || closest_channel<0) closest_channel=dist;
        }
        double furthest_same_pocket = 0; //distance to furthest MC point that is part of this pore
        vector<double> vector_of_distances_to_sphere_centroid; //keep track of the distance of each to-be-blocked point to the most dense point - it will save us calculating this again later when we check which points were blocked
        for(int j=0; j<num_selected_points; j++) { //loop over points to be blocked, and find furthest
          Point p = selected_points.at(j);
          double dist = atmnet->calcDistanceABC(p[0], p[1], p[2], most_dense[0], most_dense[1], most_dense[2]);
          if(dist>furthest_same_pocket || furthest_same_pocket<0) furthest_same_pocket=dist;
          vector_of_distances_to_sphere_centroid.push_back(dist);
        }

        //4) calculate the largest acceptable sphere at this most dense position
        double radius = 0;
        if(closest_channel<0) radius = furthest_same_pocket + probeRad + sphere_radius_overshoot; //no channels - we can safely have a single sphere, inflated by probeRad to cover the void space here
        else if(furthest_same_pocket<closest_channel) radius = min(furthest_same_pocket + probeRad + sphere_radius_overshoot, closest_channel - (probeRad + sphere_radius_overshoot)); //there are channels, but none are closer than the furthest point in this pocket, so we can extend the sphere to cover everything; try and inflate the sphere by probeRad to ensure the void space is blocked, but not if this would cause anything within probeRad of an accessible MC point to be blocked
        else radius = max(sphere_radius_overshoot, closest_channel - (probeRad + sphere_radius_overshoot)); //else there are channels closer than the furthest pocket MC point, and we have to accommodate them

        //5) create the sphere and remove to-be-blocked points that were blocked by it
        SPHERE new_sphere;
        new_sphere.x = most_dense[0];
        new_sphere.y = most_dense[1];
        new_sphere.z = most_dense[2];
        new_sphere.r = radius;
        spheres_vector.push_back(new_sphere);
        for(int j=num_selected_points-1; j>=0; j--) { //loop over points to be blocked and eliminate blocked ones - we iterate backwards so that the vector of distances is still relevant
          if(vector_of_distances_to_sphere_centroid.at(j)<radius) { //blocked
				    Point swap = selected_points.at(j);
				    selected_points.at(j) = selected_points.at(num_selected_points-1);
				    selected_points.at(num_selected_points-1) = swap;
				    selected_points.pop_back();
				    num_selected_points--; //update counter iteratively
          }
        }
      } //6) end while loop - we have totally blocked this pocket - safe to move on to next pore
    } //7) end if is_pocket - we either ran the blocking part or bypassed this pore
  } //8) end for loop - we have analysed all pores and blocked where necessary
  int num_spheres = spheres_vector.size();
  //9) write to output (directed to .block file)
  if(debug) printf("DEBUG: %d blocking spheres were created, Cartesian positions in xyz format follow, if any\n%d\nXYZ FORMAT CARTESIAN POCKETS FOR VISUALIZATION\n", num_spheres, num_spheres);
  output << num_spheres << "\n";
  for(int i=0; i<num_spheres; i++) {
    SPHERE s = spheres_vector.at(i);
    if(debug) {
      Point f; f[0] = s.x; f[1] = s.y; f[2] = s.z;
      Point c = atmnet->abc_to_xyz(f);
      printf("X %.3f %.3f %.3f %.3f\n", c[0], c[1], c[2], s.r);
    }
    output << s.x << " " << s.y << " " << s.z << " " << s.r << "\n";
  }
}

int get_most_dense_index(ATOM_NETWORK *atmnet, vector<Point> *points_vector) { //return the most dense point in a vector - if the vector is sufficiently large, only sample a subset of the points
	int most_dense_index = -1, count = 0;
	double average_dist = 0, max_density = -1;

//  int MAX_POCKET_SIZE_FOR_DENSITY_CALCULATION = 10000;
  int MAX_SAMPLES = 1000;
  int num_points = points_vector->size();
	if(num_points<1) {
		printf("ERROR: get_most_dense_index called on a vector with %d entries\n", num_points);
		exit(EXIT_FAILURE);
//} else if(num_points>=MAX_POCKET_SIZE_FOR_DENSITY_CALCULATION) { //too many points, so use sampling
//  int max_samples = MAX_POCKET_SIZE_FOR_DENSITY_CALCULATION/2;
//		double step = ((double)(num_points))/((double)(max_samples)); //step will always be at least 2
  } else {
    int max_samples = min(MAX_SAMPLES, num_points);
		double step = ((double)(num_points))/((double)(max_samples)); //step will always be at least 1
		vector<int> list;
		for(int i=0; i<max_samples; i++) {
			list.push_back((int)(((double)i)*step));
		}
		for(int i=0; i<max_samples; i++) {
			Point c1 = points_vector->at(list[i]);
			for(int j=i+1; j<max_samples; j++) {
				Point c2 = points_vector->at(list[j]);
        double dist = atmnet->calcDistanceABC(c1[0], c1[1], c1[2], c2[0], c2[1], c2[2]);
				average_dist+=dist;
				count++;
			}
		}
		average_dist/=((double)count);
		for(int i=0; i<max_samples; i++) {
			Point c1 = points_vector->at(list[i]);
			double density = 0;
			for(int j=i+1; j<max_samples; j++) {
				Point c2 = points_vector->at(list[j]);
        double dist = atmnet->calcDistanceABC(c1[0], c1[1], c1[2], c2[0], c2[1], c2[2]);
				density += exp(-1*dist*dist/(average_dist*average_dist));
			}
			if(density>max_density || max_density<0) {
				max_density = density;
				most_dense_index = i;
			}
		}
	}
/*
  else { //few enough points not to require sampling
		for(int i=0; i<num_points; i++) {
			Point c1 = points_vector->at(i);
			for(int j=i+1; j<num_points; j++) {
				Point c2 = points_vector->at(j);
        double dist = atmnet->calcDistanceABC(c1[0], c1[1], c1[2], c2[0], c2[1], c2[2]);
				average_dist+=dist;
				count++;
			}
		}
		average_dist/=((double)count);
		for(int i=0; i<num_points; i++) {
			Point c1 = points_vector->at(i);
			double density = 0;
			for(int j=i+1; j<num_points; j++) {
				Point c2 = points_vector->at(j);
        double dist = atmnet->calcDistanceABC(c1[0], c1[1], c1[2], c2[0], c2[1], c2[2]);
				density += exp(-1*dist*dist/(average_dist*average_dist));
			}
			if(density>max_density || max_density<0) {
				max_density = density;
				most_dense_index = i;
			}
		}
	}
*/
	return most_dense_index;
}










