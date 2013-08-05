/* Stores information about the basic data types that constitute a unit cell,
 * such as atoms (ATOM) and networks of atoms (ATOM_NETWORK). Also stores 
 * information about the underlying Voronoi network constituents such as nodes 
 * (VOR_NODE), edges (VOR_EDGE), and the network itself (VORONOI_NETWORK).
 */

#ifndef NETWORKSTORAGE_H
#define NETWORKSTORAGE_H

#include <iostream>
#include <string>
#include <vector>
#include "mindist.h"
#include "geometry.h"

/** Returns the determinant of the provided 3x3 matrix of doubles. */
double calcDeterminant(double matrix[3][3]);

/** Returns the determinant of the provided 3x3 matrix of integers. */
int calcDeterminant(int matrix[3][3]);

/** Store the result of inverting the provided matrix into the new
    matrix. */
void invertMatrix(double matrix [3][3], double newMatrix[3][3]);

/* Class TRIPLET and XYZ are moved to geometry.h and geometry.cc */



/** Data structure used to store the information about one atom,
    including its (x,y,z) location, its location relative to the unit
    cell axes (a_coord,b_coord,c_coord), its radius and its type. */
class ATOM {
public:
  double x,y,z;
  double a_coord, b_coord, c_coord;
  double radius;
  double mass;
  std::string type;
  int specialID;

  //added 2012-11-19 to store charge information
  double charge;
  ATOM(); //constructor sets the value to zero

  /** Print the information about this atom to the provided output stream. 
   *  Default is standard output stream*/
  void print(std::ostream &out = std::cout);

};

/** Data structure used to store information about a unit
    cell. Information includes the lengths of its three sides
    (a,b and c), its three angles (alpha, beta and gamma),
    its vectors (v_a, v_b and v_c), its name as well as all of
    its atoms. */
class ATOM_NETWORK {
 public:

  ATOM_NETWORK(); //simple constructor that sets 'false' to highAccuracyFlag

  bool highAccuracyFlag; // flag indicating if a atom network was altered (to enable high accuracy voronoi networks)
  double a, b, c;
  double alpha, beta, gamma;
  XYZ v_a, v_b, v_c;
  double ucVectors [3][3];
  double invUCVectors[3][3];
  int numAtoms;
  std::vector <ATOM> atoms ;
  std::vector <int> IDmapping;
  std::string name;
  MIN_PER_DISTANCE distanceCalculator;

  /** Print the information about this network of atoms to the
      provided output stream, including the information about each
      atom in the network. Default is standard output stream*/
  void print(std::ostream& out = std::cout);

  /** Copy the data contained in this ATOM_NETWORK to a new network using
      the provided pointer.
   */
  void copy(ATOM_NETWORK *newNet);

  /** Calculate the unit cell vectors based on the provided values
      of its side lengths and angles. v_a corresponds to the cartesian
      x-axis. */
  void initialize();

  /* Store the initialized unit cell vectors in matrix form. */
  void initMatrices();

  MIN_PER_DISTANCE getDistCalc() const;
  
    /** Determines whether a specific supercell size satisfies the non-overlapping 
     *  sphere requirement */
    //int check_sphere_overlap(int num_a, int num_b, int num_c, double diam);

    /** Identifies tetrahedra of the given atom type, calculates their tetrahedrality
     *  and returns as a vector of doubles */
    //std::vector<double> find_tetrahedra(std::string element);

    /** Determine the smallest supercell dimensions such that a sphere of a given 
     *  diameter does not overlap with itself across the periodic boundary */
    TRIPLET getSmallestSupercell(double sphere_diam);

    /** Convert coordinates relative to the unit cell (a,b,c) to absolute cartesian
    coordinates (x,y,z). */
    const Point abc_to_xyz(double a, double b, double c) const;

    /** Convert coordinates relative to the unit cell (a,b,c) to absolute cartesian
    coordinates (x,y,z). */
    const Point abc_to_xyz (Point abcPt) const;

    /** Convert coordinates relative to the unit cell (a,b,c) to absolute cartesian
    coordinates (x,y,z). */
    const Point abc_to_xyz (const XYZ& temp) const;
    
    /** Convert coordinates relative to the cartesian axes to those
    relative to the unit cell vectors. Returns a Point instance of the
    form (a_coord, b_coord, c_coord). */
    const Point xyz_to_abc(double x, double y, double z) const; 

    /** Convert coordinates relative to the cartesian axes to those
    relative to the unit cell vectors. Returns a Point instance of the
    form (a_coord, b_coord, c_coord). */
    const Point xyz_to_abc (Point xyzPt) const; 

    /** Convert coordinates relative to the cartesian axes to those
    relative to the unit cell vectors. Returns an Point instance of the
    form (a_coord, b_coord, c_coord). */
    const Point xyz_to_abc (const XYZ& temp) const;
    
    /* Calculates the minimum distance between the two points whose coordinates are relative to the unit cell vectors. */ 
    double calcDistanceABC(double a1, double b1, double c1, double a2, double b2, double c2) const;

    /** Calculates the minimum distance between the point (x,y,z) and the other point whose coordinates
     *  are relative to the unit cell vectors. */
    double calcDistanceXYZABC(double x1, double y1, double z1, double a2, double b2, double c2);

    /* Calculates the minimum distance between the point (x,y,z) and the provided atom.  */
    double calcDistance(double x, double y, double z, ATOM *atm);

    /* Calculates the minimum distance between the two provided atoms.  */
    double calcDistance(const ATOM& atm1, const ATOM& atm2) const;

    /* Calculates the minimum distance between the points (x1,y1,z1) and (x2,y2,z2).  */
    double calcDistanceXYZ(double x1, double y1, double z1, double x2, double y2, double z2);
    
    /* Rich edit: for static point (x1,y1,z1), returns the closest periodic image of point (x2,y2,z2).  */
    const XYZ getClosestPointInABC(double x1, double y1, double z1, double x2, double y2, double z2);

    /** Modify the provided (x,y,z) Point so that its coordinates reflect unit cell translations 
     *  by the provided amounts along each unit cell axis. */ 
    void translatePoint(Point *origPoint, double da, double db, double dc);

    /** Shifts the provided Point whose coordinates are relative to the unit cell vectors
     *  such that it lies within the unitcell */
    const Point shiftABCInUC(Point abcCoords);

    /** Shifts the provided Point whose coordinates are relative to the x,y,z vectors
     *  such that it lies within the unitcell */
    const Point shiftXYZInUC(Point xyzCoords);

    /** Shift the coordinates of the provided Point using the unit cell vectors until the
     *  Euclidean distance between the Point and (x,y,z) is minimal. Returns the resulting point.
     */
    const Point minimizePointDistance(Point origPoint, double dx, double dy, double dz);

    /* Identifies tetrahedra of the given atom type, calculates their tetrahedrality
     * and returns as a vector of doubles */
    std::vector<double> find_tetrahedra(std::string element);

    /* Returns Tetrahedrality index for a tetrahedra defined by four atoms */
    double CalculateTetrahedrality4Atoms(const ATOM& atm1, const ATOM& atm2, const ATOM& atm3, const ATOM& atm4) const;

};

/* Marked for deletion
/** Convert coordinates relative to the unit cell (a,b,c) to absolute cartesian
    coordinates (x,y,z). */
//Point abc_to_xyz(double a, double b, double c, ATOM_NETWORK *cell);

/** Convert coordinates relative to the unit cell (a,b,c) to absolute cartesian
    coordinates (x,y,z). */
//Point abc_to_xyz (Point abcPt, ATOM_NETWORK *cell);

/** Convert coordinates relative to the unit cell (a,b,c) to absolute cartesian
    coordinates (x,y,z). */
//Point abc_to_xyz (XYZ temp,ATOM_NETWORK *cell);

/** Convert coordinates relative to the cartesian axes to those
    relative to the unit cell vectors. Returns an Point instance of the
    form (a_coord, b_coord, c_coord). */
//Point xyz_to_abc(double x, double y, double z, ATOM_NETWORK *cell);

/** Convert coordinates relative to the cartesian axes to those
    relative to the unit cell vectors. Returns an Point instance of the
    form (a_coord, b_coord, c_coord). */
//Point xyz_to_abc (Point xyzPt, ATOM_NETWORK *cell);

/** Convert coordinates relative to the cartesian axes to those
    relative to the unit cell vectors. Returns an Point instance of the
    form (a_coord, b_coord, c_coord). */
//Point xyz_to_abc (XYZ temp, ATOM_NETWORK *cell);

/* End: Marked for deletion */


/* Determine the smallest supercell dimensions such that a sphere of a given 
   diameter does not overlap with itself across the periodic boundary */
/* Marked for deletion
   TRIPLET getSmallestSupercell(double sphere_diam, ATOM_NETWORK *atmnet);
 * End: Marked for deletion
 */

/** Data structure that stores information about an edge in a Voronoi network. */
class VOR_EDGE{
public:
  int from, to;
  double rad_moving_sphere;
  int delta_uc_x, delta_uc_y, delta_uc_z;
  double length;
  VOR_EDGE();
  VOR_EDGE(const VOR_EDGE&);
  VOR_EDGE(int myFrom, int myTo, double rad, int dx, int dy, int dz, double len);
};

/** Data structure that stores information about a node in a Voronoi network. */
class VOR_NODE{
public:
  VOR_NODE();
  VOR_NODE(double myX, double myY, double myZ, double rad, std::vector<int> ids);
  double x,y,z;
  double rad_stat_sphere;
  std::vector<int> atomIDs;
  bool active; // flag tells if a node is taken into consideration when analyzing the framework
               // e.g. this typically will set to false after pruning the network using a probe radius
};

/** Data structure that stores the nodes and edges that comprise a 
    Voronoi network. */
class VORONOI_NETWORK{
public:
    XYZ v_a, v_b, v_c;      // Unit cell vectors
    std::vector <VOR_NODE> nodes;
    std::vector <VOR_EDGE> edges;
    
    VORONOI_NETWORK();
    VORONOI_NETWORK(XYZ&, XYZ&, XYZ&, std::vector<VOR_NODE>&, std::vector<VOR_EDGE>&);

    /* Copy constructor for VORONIO_NETWORK class. */
    VORONOI_NETWORK(const VORONOI_NETWORK &);
    
    /** Copy the data contained in this VORONOI_NETWORK to a new network
     *  using the provided pointer. Deprecated, use the copy constructor.
     */
    void copy(VORONOI_NETWORK *newNet);

    /** Returns a new VORONOI_NETWORK where all the edges between nodes that are not 
     *  contained in nodeIDs are removed. However, the nodes are not removed */
    const VORONOI_NETWORK filterEdges(std::vector<int> nodeIDs);

/* Marked for deletion
    static void filterVornetEdges(std::vector<int> nodeIDs, 
                                  VORONOI_NETWORK *oldNet, 
                                  VORONOI_NETWORK *newNet);
* End: Marked for deletion
*/
    /** Copies all edges and nodes within the provided VORONOI_NETWORK
     *  to a new network iff a sphere with the specified radius can pass.*/
    //const VORONOI_NETWORK filter(const double& minRadius);

    /** Returns a copy of the VORNOI_NETWORK instance,
     *  but removes the edges that do not allow a sphere
     *  with the provided radius to pass. */
    const VORONOI_NETWORK prune(const double& minRadius);

    /** Stores a copy of the original VORNOI_NETWORK into the other provided
     *  VORONOI_NETWORK but removes the edges that are connected to specified nodes
     */
    //const VORONOI_NETWORK prunefromEdgeList(const std::vector<int>& ids);

};

/** Copies all edges and nodes within the provided VORONOI_NETWORK
    to a new network iff a sphere with the specified radius can pass.*/
void filterVoronoiNetwork(VORONOI_NETWORK *vornet, 
                          VORONOI_NETWORK *newVornet, 
                          double minRadius);


  
/** Stores a copy of the original VORNOI_NETWORK into the other provided
    VORONOI_NETWORK but removes the edges that do not allow a sphere
    with the provided radius to pass. */
/* Marked for deletion
   void pruneVoronoiNetwork(VORONOI_NETWORK *vornet, 
                         VORONOI_NETWORK *newVornet, 
                         double minRadius);
 * End: Mark for deletion
 */

/** Stores a copy of the original VORNOI_NETWORK into the other provided
 *  *  VORONOI_NETWORK but removes the edges that are connected to specified nodes
 *   */
void pruneVoronoiNetworkfromEdgeList(VORONOI_NETWORK *vornet, 
                                     VORONOI_NETWORK *newVornet, 
                                     std::vector <int> ids);

/* Attempt to substitute every other Si atom with an Al atom. ATOM_NETWORK may 
 * only consist of Si and O atoms, where each Si atom must be bonded to exactly
 * 4 oxygen atoms and each oxygen atom must be bonded to exactly 2 Si atoms. 
 * Returns true iff the substitution was successful and stores the number of 
 * substitutions using the provided reference. The provided boolean specifies 
 * whether the seeded silicon atom is substituted or not.
 * Since only 2 configurations are possible if the structure is consistent, 
 * changing this parameter enables generation of all configurations. */
bool substituteAtoms(ATOM_NETWORK *origNet, ATOM_NETWORK *newNet, 
                     bool substituteSeed, int *numSubstitutions, bool radial);

/* Attempt to substitute the specified fraction of Si atom with an Al atom. 
 * ATOM_NETWORK may only consist of Si and O atoms, where each Si atom must be 
 * bonded to exactly 4 oxygen atoms and each oxygen atom must be bonded to 
 * exactly 2 Si atoms. Returns true iff the substitution was successful and 
 * stores the number of substitutions and the  fraction of substitutions using 
 * the provided references. The provided boolean specifies whether the seeded 
 * silicon atom is substituted or not when generating the initial 50/50 
 * configuration.  The function works  by first substituting every other Si 
 * atom and then reverting some of the substituted atoms back to Si. The 
 * provided random number generator seed is used to choose which atoms to 
 * switch back.*/
bool fracSubstituteAtoms(ATOM_NETWORK *origNet, ATOM_NETWORK *newNet, 
                         bool substituteSeed, double frac, int randSeed, 
                         int *numSubstitutions, double *subFrac, bool radial);

/* Attempt to substitute the specified fraction of Si atom with an Al atom. 
 * ATOM_NETWORK may only consist of Si and O atoms, where each Si atom must be 
 * bonded to exactly 4 oxygen atoms and each oxygen atom must be bonded to 
 * exactly 2 Si atoms. 
 * This is Maciek's version that does not require initial 50/50 distribution */
bool fracSubstituteAtoms_Maciek(ATOM_NETWORK &origNet, ATOM_NETWORK &newNet, 
                                bool substituteSeed, double frac, int randSeed, 
                                int &numSubstitutions, double &subFrac, bool radial);

/* Returns the integer nearest to the provided double.*/
int nearestInt(double num);

/* Determines whether a specific supercell size satisfies the non-overlapping 
 * sphere requirement */
int check_sphere_overlap(int num_a, int num_b, int num_c, double diam, 
                         ATOM_NETWORK *atmnet);

/* Returns the id of the VOR_NODE in the provided VORONOI_NETWORK whose coordinates
 * match those of the provided Point.
 */
int getNodeID(Point pt, ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet);


#endif


