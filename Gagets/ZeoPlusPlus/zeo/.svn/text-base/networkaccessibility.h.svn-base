#ifndef ACCESSIBILITY_H
#define ACCESSIBILITY_H

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <voro++.hh>

#include "networkstorage.h"
#include "geometry.h"
#include "voronoicell.h"
#include "channel.h"

using namespace std;
using namespace voro;

/* Accessibility class handels determination of accessibility of the void space and MC-sampled points */

/* This class also handles high-accuracy calcluations where large atoms are replaced by clusters of small atoms
 when settign up this class with setupAndFindChannels() function, the first atom net provided is the high accuracy
 one based on clusters, the second is the original atom network. If the high accuracy flag is not set, both
 atom network pointers should point to the original net */

class AccessibilityClass{
    
public:
    
    ATOM_NETWORK orgAtomNet;
    ATOM_NETWORK analyzedAtomNet; // this to store the network to be analyzed (either original or the high accuracy one)
    
    bool highAccuracyFlag;
    
    VORONOI_NETWORK vornet;
    vector<BASIC_VCELL> vorcells;
    vector<VOR_CELL> advCells;
    
    vector<PORE> pores;
    //  vector<CHANNEL> channels;
    //  vector<POCKET> pockets;
    int n_channels, n_pockets;
    vector<bool> accessInfo; // flags stating if nodes are accessible
    vector<int> channelMapping; // maps node IDs to channel IDs
    vector<int> pocketMapping;  // maps node IDs to pocket IDs
    
    double r_probe;
    
    container_periodic_poly *new_rad_con;
    
    double tempMinDist; // temporary variable to store min. dis. in accessibility calcluations
    Point tempPoint; // this array stores coordiantes of the last sampled point
    int tempMinDistAtomID; // ID of the closest atom to the last investigated point
    int tempNodeID;  // ID of the voronoi node used to determine accessibility
    
    vector< pair<int, Point> > resampledInfo; // List of resampled points and the id of the Voronoi cell to which they belong
    int resampleCount;
    bool needToResampleFlag; // this flag is needed after accessibility functions have been moved to a separate class
    
public:
    
    
    AccessibilityClass(){needToResampleFlag = false; resampleCount = 0;};
    
    void setupAndFindChannels(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double r_probe_chan, double r_probe_sampl);
    
    /* checks if the provided point is in accessible volume */
    pair <bool,bool> isVPointInsideAtomAndNotAccessible(Point samplingPoint);
    
    /* checks if thr provided point on a surface of atomID is accessible */
    pair <bool,bool> isSPointInsideAtomAndNotAccessible(Point samplingPoint, int atomID);
    
    bool isVPointAccessible(Point samplingPoint);
    
    /* Checks is provide point is inside any atom and if it is accessible */
    /* the functions returns value = true if atoms overlaps/is not accessible) */
    /* atomID is an atom in the original network from where the sampling point is
     this is to check if the point is inside atom in SA calculaton */
    pair <bool,bool> isPointInsideAtomAndNotAccessible(Point samplingPoint, int atomID);
    
    /* Returns the last calculated minDist (calculated in isPointInsideAtomAndNotAccessible() )  */
    double lastMinDist(){return tempMinDist;};
    /* Returns channel or pocket ID for the last calcluated point */
    pair <int,int> lastChannelOrPocket()
    {
        if(channelMapping[tempNodeID]<0 && pocketMapping[tempNodeID]<0){
            std::cerr << "CoP_error: cannot determine point accessiblity. Consider running at higher accuracy (-ha flag)." << endl;
            //abort();
        }
        return pair <int,int> (channelMapping[tempNodeID],pocketMapping[tempNodeID]);
    };
    
    /* Returns information if a point need to be resampled due to node not found in isPointInsideAtomAndNotAccessible() */
    bool needToResample()
    {
        if(needToResampleFlag) {
            needToResampleFlag = false;
            return true;
        }
        else return false;
    }
    
    /* Return number of resampled ponits */
    int getResampleCount(){return resampleCount;};
    
    /* remove Voronoi nodes that are not used in analysis */
    void removeOverlappedNodes();
    
    /* deconstruct class */
    void deconstruct(){delete new_rad_con;};
};


#endif

