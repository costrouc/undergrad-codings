/* Stores information about the basic data types that constitute a unit cell,
 * such as atoms (ATOM) and networks of atoms (ATOM_NETWORK). Also stores 
 * information about the underlying Voronoi network constituents such as nodes 
 * (VOR_NODE), edges (VOR_EDGE), and the network itself (VORONOI_NETWORK).
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <map>
#include <algorithm>
#include <float.h>
#include "networkstorage.h"
#include "networkinfo.h"
#include "geometry.h"
#include "zeo_consts.h"
using namespace std;

ATOM_NETWORK::ATOM_NETWORK() {highAccuracyFlag = false;}

/** Returns the determinant of the provided 3x3 matrix of doubles. */
double calcDeterminant(double matrix[3][3]){
  return   
      matrix[0][0]*(matrix[2][2]*matrix[1][1] - matrix[2][1]*matrix[1][2])
    - matrix[1][0]*(matrix[2][2]*matrix[0][1] - matrix[2][1]*matrix[0][2])
    + matrix[2][0]*(matrix[1][2]*matrix[0][1] - matrix[1][1]*matrix[0][2]);
}

/** Returns the determinant of the provided 3x3 matrix of integers. */
int calcDeterminant(int matrix[3][3]){
  return   
    matrix[0][0]*(matrix[2][2]*matrix[1][1] - matrix[2][1]*matrix[1][2])
    - matrix[1][0]*(matrix[2][2]*matrix[0][1] - matrix[2][1]*matrix[0][2])
    + matrix[2][0]*(matrix[1][2]*matrix[0][1] - matrix[1][1]*matrix[0][2]);
}

/** Store the result of inverting the provided matrix into the new
    matrix. */
void invertMatrix(double matrix [3][3], double newMatrix[3][3]){
  double determinant = calcDeterminant(matrix);
  if(determinant == 0){
    cerr << "Determinant of provided matrix is 0. Matrix is not invertible."
	 << "\n" << "Exiting ..." << "\n";
    exit(1);
  }
    
  double invDet = 1/determinant;
  newMatrix[0][0] = invDet*   (matrix[2][2]*matrix[1][1]-matrix[2][1]*matrix[1][2]);
  newMatrix[0][1] = invDet*-1*(matrix[2][2]*matrix[0][1]-matrix[2][1]*matrix[0][2]);
  newMatrix[0][2] = invDet*   (matrix[1][2]*matrix[0][1]-matrix[1][1]*matrix[0][2]);
  newMatrix[1][0] = invDet*-1*(matrix[2][2]*matrix[1][0]-matrix[2][0]*matrix[1][2]);
  newMatrix[1][1] = invDet*   (matrix[2][2]*matrix[0][0]-matrix[2][0]*matrix[0][2]);
  newMatrix[1][2] = invDet*-1*(matrix[1][2]*matrix[0][0]-matrix[1][0]*matrix[0][2]);
  newMatrix[2][0] = invDet*   (matrix[2][1]*matrix[1][0]-matrix[2][0]*matrix[1][1]);
  newMatrix[2][1] = invDet*-1*(matrix[2][1]*matrix[0][0]-matrix[2][0]*matrix[0][1]);
  newMatrix[2][2] = invDet*   (matrix[1][1]*matrix[0][0]-matrix[1][0]*matrix[0][1]);
}



ATOM::ATOM() {
  charge = 0; //because this is not always provided, a default value of zero is set
}

/** Print the information about this atom to the provided output stream. */
void ATOM::print(ostream &out){
  out << "   type:" << type << "   x:" << x << "   y:" << y << "   z:"  << z  
      << "   a:" << a_coord << "   b:" << b_coord << "   c:" << c_coord
      << "   radius:" << radius << "\n";
}


/** Print the information about this network of atoms to the
 *  provided output stream, including the information about each
 *  atom in the network. */
void ATOM_NETWORK::print(ostream &out){
  out << "Name: " << name << "\n"
      << "A: " << a <<"     B: " << b << "     C: " << c << "\n"
      << "Alpha: " << alpha <<"     Beta: " << beta << "     Gamma: " << gamma 
      << "\n";
  out << "v_a: "; v_a.print();
  out << "v_b: "; v_b.print();
  out << "v_c: "; v_c.print();
  out << "Number of atoms: " << numAtoms << "\n";
  for(int i = 0; i<numAtoms; i++)
    atoms.at(i).print();
}


/** Copy the data contained in this ATOM_NETWORK to a new network using
    the provided pointer.
*/
void ATOM_NETWORK::copy(ATOM_NETWORK *newNet){
  newNet->a = a; newNet->b = b; newNet->c = c;
  newNet->alpha = alpha; newNet->beta = beta; newNet->gamma = gamma;
  newNet->v_a = v_a; newNet->v_b = v_b; newNet->v_c = v_c;
  newNet->numAtoms = numAtoms; newNet->name = name;
  newNet->atoms.clear(); newNet->atoms = atoms;
  newNet->IDmapping.clear(); newNet->IDmapping = IDmapping;
  newNet->initialize();
}


/** Calculate the unit cell vectors based on the provided values
    of its side lengths and angles. v_a corresponds to the cartesian
    x-axis. 
*/
void ATOM_NETWORK::initialize(){
  double tempd, talpha, tbeta, tgamma;
  talpha = 2*PI/360.0*alpha;
  tbeta  = 2*PI/360.0*beta;
  tgamma = 2*PI/360.0*gamma;
  tempd=(cos(talpha)-cos(tgamma)*cos(tbeta))/sin(tgamma);
  
  v_a.x=a;
  v_a.y=0.0;
  v_a.z=0.0;
  
  v_b.x=b*cos(tgamma); 
  if(fabs(v_b.x) < TOLERANCE)
    v_b.x=0;
  
  v_b.y=b*sin(tgamma);
  v_b.z=0.0;
  
  v_c.x=c*cos(tbeta);
  if(fabs(v_c.x) < TOLERANCE)
    v_c.x = 0;
  v_c.y=c*tempd;
  if(fabs(v_c.y) < TOLERANCE)
    v_c.y = 0;
  v_c.z=c*sqrt(1.0-SQR(cos(tbeta))-SQR(tempd));
  initMatrices();
  distanceCalculator = MIN_PER_DISTANCE(v_a.x, v_b.x, v_b.y, v_c.x, v_c.y, v_c.z);
}

/* Store the initialized unit cell vectors in matrix form. */
void ATOM_NETWORK::initMatrices(){
  ucVectors[0][0] = v_a.x; ucVectors[1][0] = v_a.y; ucVectors[2][0] = v_a.z;
  ucVectors[0][1] = v_b.x; ucVectors[1][1] = v_b.y; ucVectors[2][1] = v_b.z;
  ucVectors[0][2] = v_c.x; ucVectors[1][2] = v_c.y; ucVectors[2][2] = v_c.z;
  invertMatrix(ucVectors,invUCVectors);
}

MIN_PER_DISTANCE ATOM_NETWORK::getDistCalc() const {
  return distanceCalculator;
}


/** Determine the smallest supercell dimensions such that a sphere of a given 
 *  diameter does not overlap with itself across the periodic boundary */
TRIPLET ATOM_NETWORK::getSmallestSupercell(double diam) {
  //1) set lower bound based on lengths of cell axes
  int na = 1+ (diam/a);
  int nb = 1+ (diam/b);
  int nc = 1+ (diam/c);

  int fewest_cells = -1; //no supercell yet successful
  TRIPLET smallest_supercell(-1,-1,-1);

  //2) search all possible supercell sizes to find the smallest one satisfying 
  // the minimum image convention for this radius
  TRIPLET lb(na, nb, nc);
  vector<TRIPLET> supercells;
  supercells.push_back(lb);
  while(supercells.size()>0) {
    //3) take the oldest candidate on the vector
    TRIPLET s = supercells.at(0); //oldest supercell candidate, s
    for(int i=0; i<supercells.size()-1; i++) 
      supercells.at(i) = supercells.at(i+1); //shift vector up
    supercells.pop_back(); //delete last, which is now duplicated at last-1
    int num_cells = s.x*s.y*s.z;
    //4) is s a potential new smallest supercell?
    if(num_cells<fewest_cells || fewest_cells<0) {
      //5) does s satisfy the minimum image convention?
      int status = check_sphere_overlap(s.x, s.y, s.z, diam, this); //For time being
      if(status==-1) {
        printf("WARNING: bad unit cell angles!\n");
        return smallest_supercell;
      } else if(status==1) { //acceptable!
        fewest_cells = num_cells;
        smallest_supercell = s;
        //printf("smallest satisfactory supercell so far: (%d %d %d) = 
                  //%d cells\n", s.a, s.b, s.c, num_cells);
      } else { //unacceptable - try larger supercells in each direction
        TRIPLET s2(s.x+1, s.y, s.z);
        TRIPLET s3(s.x, s.y+1, s.z);
        TRIPLET s4(s.x, s.y, s.z+1);
        supercells.push_back(s2);
        supercells.push_back(s3);
        supercells.push_back(s4);
      }
    }
  }
  return smallest_supercell;
}


// abc_to_xyz and xyz_to_abc functions which need ATOM_NETWORK are moved from geometry
/** Convert coordinates relative to the unit cell (a,b,c) to absolute cartesian
    coordinates (x,y,z). */
const Point ATOM_NETWORK::abc_to_xyz(double a, double b, double c) const {
  //Point xyzCoords;
  /*
  xyzCoords.x = a*cell->v_a.x+b*cell->v_b.x+c*cell->v_c.x;
  xyzCoords.y = a*cell->v_a.y+b*cell->v_b.y+c*cell->v_c.y;
  xyzCoords.z = a*cell->v_a.z+b*cell->v_b.z+c*cell->v_c.z;
  if(fabs(xyzCoords.x) < 0.00001) xyzCoords.x = 0;
  if(fabs(xyzCoords.y) < 0.00001) xyzCoords.y = 0;
  if(fabs(xyzCoords.z) < 0.00001) xyzCoords.z = 0;
  */

  // Use only non-zero elements in computation
  double xt = a*v_a.x+b*v_b.x+c*v_c.x;
  double yt = b*v_b.y+c*v_c.y;
  double zt = c*v_c.z;
  return Point(xt, yt, zt);
}

/** Convert coordinates relative to the unit cell (a,b,c) to absolute cartesian
    coordinates (x,y,z). */
const Point ATOM_NETWORK::abc_to_xyz(Point abcPt) const {
  return abc_to_xyz(abcPt[0], abcPt[1], abcPt[2]);
}

const Point ATOM_NETWORK::abc_to_xyz (const XYZ& temp) const {
  return abc_to_xyz(temp.x,temp.y,temp.z);
}

/** Convert coordinates relative to the cartesian axes to those
    relative to the unit cell vectors. Returns an Point instance of the
    form (a_coord, b_coord, c_coord). */
const Point ATOM_NETWORK::xyz_to_abc(double xi, double yi, double zi) const{
  //Point abcCoords;
  /*
  abcCoords.x = x*cell->invUCVectors[0][0]+y*cell->invUCVectors[0][1]+z*cell->invUCVectors[0][2];
  abcCoords.y = x*cell->invUCVectors[1][0]+y*cell->invUCVectors[1][1]+z*cell->invUCVectors[1][2];
  abcCoords.z = x*cell->invUCVectors[2][0]+y*cell->invUCVectors[2][1]+z*cell->invUCVectors[2][2];
  if(fabs(abcCoords.x) < 0.0001) abcCoords.x = 0;
  if(fabs(abcCoords.y) < 0.0001) abcCoords.y = 0;
  if(fabs(abcCoords.z) < 0.0001) abcCoords.z = 0;
  */

  // Use only non-zero elements in computation
  double xt = xi*invUCVectors[0][0]+yi*invUCVectors[0][1]+zi*invUCVectors[0][2];
  double yt = yi*invUCVectors[1][1]+zi*invUCVectors[1][2];
  double zt = zi*invUCVectors[2][2];
  return Point(xt, yt, zt);
}

/** Convert coordinates relative to the cartesian axes to those
    relative to the unit cell vectors. Returns an Point instance of the
    form (a_coord, b_coord, c_coord). */
const Point ATOM_NETWORK::xyz_to_abc (Point xyzPt) const {
  return xyz_to_abc(xyzPt[0], xyzPt[1], xyzPt[2]);
}

const Point ATOM_NETWORK::xyz_to_abc(const XYZ& temp) const {
  return xyz_to_abc(temp.x,temp.y,temp.z);
}

/** Calculates the minimum distance between the two points whose coordinates are relative to the unit cell vectors.    */ 
double ATOM_NETWORK::calcDistanceABC(double a1, double b1, double c1, double a2, double b2, double c2) const{
  return getDistCalc().minimum_periodic_distance(a1, b1, c1, a2, b2, c2);
}

/** Calculates the minimum distance between the point (x,y,z) and the other point whose coordinates
 *  are relative to the unit cell vectors. */
double ATOM_NETWORK::calcDistanceXYZABC(double x1, double y1, double z1, double a2, double b2, double c2){
  Point abcCoord = xyz_to_abc(x1, y1, z1);
  return calcDistanceABC(abcCoord[0], abcCoord[1], abcCoord[2], a2, b2, c2);
}

/** Calculates the minimum distance between the point (x,y,z) and the provided atom.  */
double ATOM_NETWORK::calcDistance(double x, double y, double z, ATOM *atm){
  return calcDistanceXYZABC(x, y, z, atm->a_coord, atm->b_coord, atm->c_coord);
}

/** Calculates the minimum distance between the two provided atoms.  */
double ATOM_NETWORK::calcDistance(const ATOM& atm1, const ATOM& atm2) const {
  return calcDistanceABC(atm1.a_coord, atm1.b_coord, atm1.c_coord, atm2.a_coord, atm2.b_coord, atm2.c_coord); 
}

/** Calculates the minimum distance between the points (x1,y1,z1) and (x2,y2,z2).  */
double ATOM_NETWORK::calcDistanceXYZ(double x1, double y1, double z1, double x2, double y2, double z2){
  Point abcCoord = xyz_to_abc(x1, y1, z1);
  return calcDistanceXYZABC(x2, y2, z2, abcCoord[0], abcCoord[1], abcCoord[2]);
}

/** Rich edit: for static point (x1,y1,z1), returns the closest periodic image of point (x2,y2,z2).*/ 
const XYZ ATOM_NETWORK::getClosestPointInABC(double x1, double y1, double z1, double x2, double y2, double z2){
  Point abcCoordStatic = xyz_to_abc(x1, y1, z1);
  Point abcCoordMobile = xyz_to_abc(x2, y2, z2);
  XYZ answer;
  getDistCalc().closest_periodic_image(abcCoordStatic[0], abcCoordStatic[1], abcCoordStatic[2], 
					       abcCoordMobile[0], abcCoordMobile[1], abcCoordMobile[2], 
					       answer.x, answer.y, answer.z);
  return answer;
}

/** Modify the provided (x,y,z) Point so that its coordinates reflect unit cell translations 
 *  by the provided amounts along each unit cell axis. */ 
void ATOM_NETWORK::translatePoint(Point *origPoint, double da, double db, double dc){
  (*origPoint)[0] = (*origPoint)[0] +  da*v_a.x + db*v_b.x + dc*v_c.x;
  (*origPoint)[1] = (*origPoint)[1] +  da*v_a.y + db*v_b.y + dc*v_c.y;
  (*origPoint)[2] = (*origPoint)[2] +  da*v_a.z + db*v_b.z + dc*v_c.z;
}

/** Shifts the provided Point whose coordinates are relative to the unit cell vectors
 *  such that it lies within the unitcell.  */
const Point ATOM_NETWORK::shiftABCInUC(Point abcCoords) {
  return Point(trans_to_origuc(abcCoords[0]), trans_to_origuc(abcCoords[1]), trans_to_origuc(abcCoords[2]));
}

/** Shifts the provided Point whose coordinates are relative to the x,y,z vectors
 *  such that it lies within the unitcell   */
const Point ATOM_NETWORK::shiftXYZInUC(Point xyzCoords) {
  Point abcCoords = shiftABCInUC(xyz_to_abc(xyzCoords));
  return abc_to_xyz(abcCoords[0], abcCoords[1], abcCoords[2]);
}

/** Shift the coordinates of the provided Point using the unit cell vectors until
 *  the Euclidean distance between the Point and (x,y,z) is minimal. Returns the resulting point.
 */
const Point ATOM_NETWORK::minimizePointDistance(Point origPoint, double dx, double dy, double dz){
  Point abc_one  = xyz_to_abc(origPoint);
  Point abc_two  = xyz_to_abc(dx, dy, dz);
  
  double minDa = DBL_MAX, minDb = DBL_MAX, minDc = DBL_MAX, best_a = DBL_MAX, best_b = DBL_MAX, best_c = DBL_MAX;
  getDistCalc().closest_periodic_image(abc_two[0], abc_two[1], abc_two[2], abc_one[0], abc_one[1], abc_one[2], 
					       minDa, minDb, minDc, best_a, best_b, best_c);
  return abc_to_xyz(best_a, best_b, best_c);
}

/* Identifies tetrahedra of the given atom type, calculates their tetrahedrality and returns as a vector of doubles */
vector<double> ATOM_NETWORK::find_tetrahedra(string element) {
  vector<double> tetras;
  double d_ij = -1, d_ik = -1, d_il = -1, d_jk = -1, d_jl = -1, d_kl = -1; //stores distances between each pair of atoms in a potential tetrahedral arrangement
  double maxDist = 5, minDist = 0.1; //four atoms are considered to be in a tetrahedral arrangement if they are each within this distance range of each other - if this is true, we then calculate the index of tetrahedral distortion
  for(int i=0; i<numAtoms; i++) {
    if(atoms.at(i).type.compare(element)==0) {
      for(int j=i+1; j<numAtoms; j++) {
        if(atoms.at(j).type.compare(element)==0) {
          d_ij = getDistCalc().minimum_periodic_distance(atoms.at(i).a_coord, atoms.at(i).b_coord, atoms.at(i).c_coord, atoms.at(j).a_coord, atoms.at(j).b_coord, atoms.at(j).c_coord);
          if(d_ij>minDist && d_ij<maxDist) {
            for(int k=j+1; k<numAtoms; k++) {
              if(atoms.at(k).type.compare(element)==0) {
                d_ik = getDistCalc().minimum_periodic_distance(atoms.at(i).a_coord, atoms.at(i).b_coord, atoms.at(i).c_coord, atoms.at(k).a_coord, atoms.at(k).b_coord, atoms.at(k).c_coord);
                if(d_ik>minDist && d_ik<maxDist) {
                  d_jk = getDistCalc().minimum_periodic_distance(atoms.at(j).a_coord, atoms.at(j).b_coord, atoms.at(j).c_coord, atoms.at(k).a_coord, atoms.at(k).b_coord, atoms.at(k).c_coord);
                  if(d_jk>minDist && d_jk<maxDist) {
                    for(int l=k+1; l<numAtoms; l++) {
                      if(atoms.at(l).type.compare(element)==0) {
                        d_il = getDistCalc().minimum_periodic_distance(atoms.at(i).a_coord, atoms.at(i).b_coord, atoms.at(i).c_coord, atoms.at(l).a_coord, atoms.at(l).b_coord, atoms.at(l).c_coord);
                        if(d_il>minDist && d_il<maxDist) {
                          d_jl = getDistCalc().minimum_periodic_distance(atoms.at(j).a_coord, atoms.at(j).b_coord, atoms.at(j).c_coord, atoms.at(l).a_coord, atoms.at(l).b_coord, atoms.at(l).c_coord);
                          if(d_jl>minDist && d_jl<maxDist) {
                            d_kl = getDistCalc().minimum_periodic_distance(atoms.at(k).a_coord, atoms.at(k).b_coord, atoms.at(k).c_coord, atoms.at(l).a_coord, atoms.at(l).b_coord, atoms.at(l).c_coord);
                            if(d_kl>minDist && d_kl<maxDist) {
                              //we have found four atoms within the tolerance
                              double t = CalculateTetrahedrality4Atoms(atoms.at(i), atoms.at(j), atoms.at(k), atoms.at(l));
                              tetras.push_back(t);
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  //sort before returning
  sort(tetras.begin(), tetras.end());
  return tetras;
}

/* Returns Tetrahedrality index for a tetrahedron defined by four atoms */
double ATOM_NETWORK::CalculateTetrahedrality4Atoms(const ATOM& atm1, const ATOM& atm2, const ATOM& atm3, const ATOM& atm4) const {

  vector <double> edges; // tetrahedra edges' lengths
  double d_mean=0.0; // mean value
  double Tindex=0.0;

  edges.push_back(calcDistance(atm1,atm2));
  edges.push_back(calcDistance(atm1,atm3));
  edges.push_back(calcDistance(atm1,atm4));
  edges.push_back(calcDistance(atm2,atm3));
  edges.push_back(calcDistance(atm2,atm4));
  edges.push_back(calcDistance(atm3,atm4));

  for(int i=0;i<6;i++) d_mean+=edges[i];

  d_mean=d_mean/6.0;

  for(int i=0;i<5;i++) {
    for(int j=i+1;j<6;j++) {
      Tindex+=((edges[i]-edges[j])*(edges[i]-edges[j]))/(15.0*d_mean*d_mean);
    }
  }

  return Tindex;
}



/* Determine the smallest supercell dimensions such that a sphere of a given 
 * diameter does not overlap with itself across the periodic boundary */
/* Marked for deletion
TRIPLET getSmallestSupercell(double diam, ATOM_NETWORK *atmnet) {
  //1) set lower bound based on lengths of cell axes
  int na = 1+ (diam/atmnet->a);
  int nb = 1+ (diam/atmnet->b);
  int nc = 1+ (diam/atmnet->c);

  int fewest_cells = -1; //no supercell yet successful
  TRIPLET smallest_supercell(-1,-1,-1);

  //2) search all possible supercell sizes to find the smallest one satisfying 
  // the minimum image convention for this radius
  TRIPLET lb(na, nb, nc);
  vector<TRIPLET> supercells;
  supercells.push_back(lb);
  while(supercells.size()>0) {
    //3) take the oldest candidate on the vector
    TRIPLET s = supercells.at(0); //oldest supercell candidate, s
    for(int i=0; i<supercells.size()-1; i++) 
      supercells.at(i) = supercells.at(i+1); //shift vector up
    supercells.pop_back(); //delete last, which is now duplicated at last-1
    int num_cells = s.x*s.y*s.z;
    //4) is s a potential new smallest supercell?
    if(num_cells<fewest_cells || fewest_cells<0) {
      //5) does s satisfy the minimum image convention?
      int status = check_sphere_overlap(s.x, s.y, s.z, diam, atmnet);
      if(status==-1) {
        printf("WARNING: bad unit cell angles!\n");
        return smallest_supercell;
      } else if(status==1) { //acceptable!
        fewest_cells = num_cells;
        smallest_supercell = s;
        //printf("smallest satisfactory supercell so far: (%d %d %d) = 
                  //%d cells\n", s.a, s.b, s.c, num_cells);
      } else { //unacceptable - try larger supercells in each direction
        TRIPLET s2(s.x+1, s.y, s.z);
        TRIPLET s3(s.x, s.y+1, s.z);
        TRIPLET s4(s.x, s.y, s.z+1);
        supercells.push_back(s2);
        supercells.push_back(s3);
        supercells.push_back(s4);
      }
    }
  }
  return smallest_supercell;
}
* End: Marked for deletion
*/

/// Default constructor
VOR_EDGE::VOR_EDGE(){}

/// Copy constructor
VOR_EDGE::VOR_EDGE(const VOR_EDGE& orig) :
    from(orig.from), to(orig.to), rad_moving_sphere(orig.rad_moving_sphere),
    delta_uc_x(orig.delta_uc_x), delta_uc_y(orig.delta_uc_y),
    delta_uc_z(orig.delta_uc_z) {}
         
VOR_EDGE::VOR_EDGE(int myFrom, int myTo, double rad, 
                   int dx, int dy, int dz, double len){
  from = myFrom; to = myTo; rad_moving_sphere = rad;
  delta_uc_x = dx;
  delta_uc_y = dy;
  delta_uc_z = dz;
  length = len;
}



VOR_NODE::VOR_NODE(){}
VOR_NODE::VOR_NODE(double myX, double myY, double myZ, 
                   double rad, vector<int> ids){
  x = myX; y = myY; z = myZ;
  rad_stat_sphere = rad;
  atomIDs = ids;
}


/// Voronoi network constructor
VORONOI_NETWORK::VORONOI_NETWORK (){}
VORONOI_NETWORK::VORONOI_NETWORK (XYZ& inp_va, XYZ& inp_vb, XYZ& inp_vc, 
                vector<VOR_NODE>& inp_nodes, vector<VOR_EDGE>& inp_edges)
                : edges(inp_edges), nodes(inp_nodes), v_a(inp_va),
                v_b(inp_vb), v_c(inp_vc) {}

/** Copy the data contained in this VORONOI_NETWORK to a new network using
    the provided pointer. Deprecated, use copy constructor.
*/
void VORONOI_NETWORK::copy(VORONOI_NETWORK *newNet){
  newNet->v_a = v_a; newNet->v_b = v_b; newNet->v_c = v_c;
  newNet->edges.clear(); newNet->edges = edges;
  newNet->nodes.clear(); newNet->nodes = nodes;
}

/// Copy constructor for VORONOI_NETWORK
VORONOI_NETWORK::VORONOI_NETWORK (const VORONOI_NETWORK& oldNet)
     : v_a(oldNet.v_a), v_b(oldNet.v_b), v_c(oldNet.v_c),
     edges(oldNet.edges), nodes(oldNet.nodes) {}



/* Marked for deletion
void VORONOI_NETWORK::filterVornetEdges(vector<int> nodeIDs, 
                                        VORONOI_NETWORK *oldNet, 
                                        VORONOI_NETWORK *newNet)
{    
    vector<bool> includeNodes = vector<bool>(oldNet->nodes.size(), false);
    for(unsigned int i = 0; i < nodeIDs.size(); i++){
      includeNodes[nodeIDs[i]] = true;
    }
    
    vector<VOR_NODE> newNodes = vector<VOR_NODE>();
    for(unsigned int i = 0; i < oldNet->nodes.size(); i++){
	newNodes.push_back(oldNet->nodes[i]);
    }
    
    vector<VOR_EDGE> newEdges = vector<VOR_EDGE>();
    
    for(unsigned int i = 0; i < oldNet->edges.size(); i++){
      VOR_EDGE edge = oldNet->edges[i];
      if(includeNodes[edge.from] && includeNodes[edge.to]){
	newEdges.push_back(edge);
      }
    }
    
    newNet->nodes = newNodes;
    newNet->edges = newEdges;
    newNet->v_a = oldNet->v_a; 
    newNet->v_b = oldNet->v_b; 
    newNet->v_c = oldNet->v_c;
}
* End: marked for deletion
*/

/* Removes all the edges between nodes that are not both contained in nodeIDs. 
   However, these nodes remain in the Voronoi network 
*/
const VORONOI_NETWORK VORONOI_NETWORK::filterEdges(vector<int> nodeIDs)
{    
    vector<bool> includeNodes = vector<bool>(nodes.size(), false);
    for(unsigned int i = 0; i < nodeIDs.size(); i++){
        includeNodes[nodeIDs[i]] = true;
    }
    
    vector<VOR_NODE> newNodes = vector<VOR_NODE>();
    for(unsigned int i = 0; i < nodes.size(); i++){
        newNodes.push_back(nodes[i]);
    }
    
    vector<VOR_EDGE> newEdges = vector<VOR_EDGE>();
    for(unsigned int i = 0; i < edges.size(); i++){
        VOR_EDGE edge = edges[i];
        if(includeNodes[edge.from] && includeNodes[edge.to]){
            newEdges.push_back(edge);
        }
    }

    return VORONOI_NETWORK(v_a, v_b, v_c, newNodes, newEdges);
}

/** Copies all edges and nodes within the provided VORONOI_NETWORK
    to a new network iff a sphere with the specified radius can pass.*/
void filterVoronoiNetwork(VORONOI_NETWORK *vornet, 
                          VORONOI_NETWORK *newVornet, 
                          double minRadius){
   map<int,int> idMappings;
  vector<VOR_NODE>::iterator nodeIter = vornet->nodes.begin();
  vector<VOR_NODE> newNodes;
  int i = 0;
  int idCount = 0;

  // Add all nodes whose radius is greater than the provided minimum to a list
  while(nodeIter != vornet->nodes.end()){
    if(nodeIter->rad_stat_sphere > minRadius){
      newNodes.push_back(*nodeIter);
      idMappings.insert(pair<int,int> (i, idCount));
      idCount++;
    }
    i++;
    nodeIter++;
  }

  // Copy nodes that met requirement
  newVornet->nodes = newNodes;

  vector<VOR_EDGE>::iterator edgeIter = vornet->edges.begin();
  vector<VOR_EDGE> newEdges;

   // Add all edges whose radius is greater than the provided minimum to a list
  while(edgeIter != vornet->edges.end()){
    if((edgeIter->rad_moving_sphere > minRadius) && (
        idMappings.find(edgeIter->from) != idMappings.end()) && (
        idMappings.find(edgeIter->to) != idMappings.end())){
      VOR_EDGE newEdge;
      newEdge.from = idMappings.find(edgeIter->from)->second; 
      newEdge.to   = idMappings.find(edgeIter->to)->second; 
      newEdge.rad_moving_sphere = edgeIter->rad_moving_sphere;
      newEdge.delta_uc_x = edgeIter->delta_uc_x; 
      newEdge.delta_uc_y = edgeIter->delta_uc_y; 
      newEdge.delta_uc_z = edgeIter->delta_uc_z; 
      newEdge.length = edgeIter->length;
      newEdges.push_back(newEdge);
    }
    edgeIter++;
  } 

   // Copy edges that met requirement
  newVornet->edges = newEdges;

  // Copy unitcell vectors to new network
  newVornet->v_a = vornet->v_a;
  newVornet->v_b = vornet->v_b;
  newVornet->v_c = vornet->v_c;
}


/** Returns a copy of the VORNOI_NETWORK instance 
    but removes the edges that do not allow a sphere
    with the provided radius to pass. */
const VORONOI_NETWORK VORONOI_NETWORK::prune(const double& minRadius)
{
    vector<VOR_EDGE> newEdges;
    // Add edges whose radius is greater than the input minimum to a list
    for(vector<VOR_EDGE>::iterator edgeIter = edges.begin();
            edgeIter != edges.end(); edgeIter++) {
        if(edgeIter->rad_moving_sphere > minRadius){
            //VOR_EDGE newEdge = *edgeIter;
            newEdges.push_back(*edgeIter);
        }
    };
   
    vector<VOR_NODE> newNodes = nodes;
    for(unsigned int i = 0; i < nodes.size(); i++)
       {
       if(nodes[i].rad_stat_sphere > minRadius) newNodes[i].active = true; else newNodes[i].active = false;
       }; 

    return VORONOI_NETWORK(v_a, v_b, v_c, newNodes, newEdges);
}


/** Stores a copy of the original VORNOI_NETWORK into the other provided
    VORONOI_NETWORK but removes the edges that do not allow a sphere
    with the provided radius to pass. */
/* Marked for deletion
void pruneVoronoiNetwork(VORONOI_NETWORK *vornet, 
                         VORONOI_NETWORK *newVornet, 
                         double minRadius){
  newVornet->nodes = vornet->nodes;

  vector<VOR_EDGE>::iterator edgeIter = vornet->edges.begin();
  vector<VOR_EDGE> newEdges;

  // Add all edges whose radius is greater than the provided minimum to a list
  while(edgeIter != vornet->edges.end()){
    if(edgeIter->rad_moving_sphere > minRadius){
      VOR_EDGE newEdge;
      newEdge.from = edgeIter->from; 
      newEdge.to   = edgeIter->to;
      newEdge.rad_moving_sphere = edgeIter->rad_moving_sphere;
      newEdge.delta_uc_x = edgeIter->delta_uc_x; 
      newEdge.delta_uc_y = edgeIter->delta_uc_y; 
      newEdge.delta_uc_z = edgeIter->delta_uc_z; 
      newEdge.length = edgeIter->length;
      newEdges.push_back(newEdge);
    }
    edgeIter++;
  } 

  // Copy edges that met requirement
  newVornet->edges = newEdges;

  // Copy unitcell vectors to new network
  newVornet->v_a = vornet->v_a;
  newVornet->v_b = vornet->v_b;
  newVornet->v_c = vornet->v_c;
}
* End: Marked for deletion
*/

/** Stores a copy of the original VORNOI_NETWORK into the other provided
 *  VORONOI_NETWORK but removes the edges that are connected to specified 
 *  nodes (specified by ID list)
 */
void pruneVoronoiNetworkfromEdgeList(VORONOI_NETWORK *vornet, 
                                     VORONOI_NETWORK *newVornet, 
                                     vector <int> ids){

  newVornet->nodes = vornet->nodes;

  vector<VOR_EDGE>::iterator edgeIter = vornet->edges.begin();
  vector<VOR_EDGE> newEdges;

 // Add all edges whose are not connected to specified nodes
  while(edgeIter != vornet->edges.end()){
    int flag=0;
    for(unsigned int i=0; i<ids.size(); i++){
	  // if edge connects to node of specified id, flag it
	  if(edgeIter->from==ids[i]||edgeIter->to==ids[i]) flag++;  
	};

    if(flag==0){ // only keep unflagged edges
         VOR_EDGE newEdge;
         newEdge.from = edgeIter->from;
         newEdge.to   = edgeIter->to;
         newEdge.rad_moving_sphere = edgeIter->rad_moving_sphere;
         newEdge.delta_uc_x = edgeIter->delta_uc_x; 
         newEdge.delta_uc_y = edgeIter->delta_uc_y;
         newEdge.delta_uc_z = edgeIter->delta_uc_z; 
         newEdge.length = edgeIter->length;
         newEdges.push_back(newEdge);
    }
    edgeIter++;
  }
  

  newVornet->edges = newEdges;
  newVornet->v_a = vornet->v_a;
  newVornet->v_b = vornet->v_b;
  newVornet->v_c = vornet->v_c;
}


/* Attempt to substitute every other Si atom with an Al atom. ATOM_NETWORK may 
 * only consist of Si and O atoms, where each Si atom must be bonded to exactly
 * 4 oxygen atoms and each oxygen atom must be bonded to exactly 2 Si atoms. 
 * Returns true iff the substitution was successful and stores the number of 
 * substitutions using the provided reference. The provided boolean specifies 
 * whether the seeded silicon atom is substituted or not.
 * Since only 2 configurations are possible if the structure is consistent, 
 * changing this parameter enables generation of all configurations. */
bool substituteAtoms(ATOM_NETWORK *origNet, ATOM_NETWORK *newNet, 
                     bool substituteSeed, int *numSubstitutions, bool radial){
  int numAtoms = origNet->numAtoms;
  double max_bond_length = 1.95;

  vector< vector<int> > bonds = vector< vector<int> >(numAtoms, vector<int>());
  for(int i = 0; i < numAtoms; i++){
    ATOM atom_one = origNet->atoms[i];
    for(int j = i + 1; j < numAtoms; j++){
      ATOM atom_two = origNet->atoms[j];
      if(origNet->calcDistance(atom_one, atom_two) < max_bond_length){
	if(atom_one.type.compare(atom_two.type) == 0){
	  cerr << "Atomic network substitution aborted because atoms of same type "
	          "are bonded to one another" << "\n" << "Occurred for type " 
	       << atom_one.type << " between atoms " << i << " and " << j << "\n";
	  return false;
	}
	bonds[i].push_back(j);
	bonds[j].push_back(i);
      }									    
    }
  }

  for(int i = 0; i < numAtoms; i++){
    if(origNet->atoms[i].type.compare("Si") == 0){
      if(bonds[i].size() != 4){
	cerr << "Atomic network substitution aborted because Si atom bonded to " 
	     << bonds[i].size() << " other atoms instead of 4" << "\n"
	     << "Occurred for atom " << i << "\n";
	return false;
      }
    }
    else if(origNet->atoms[i].type.compare("O") == 0){
      if(bonds[i].size() != 2){
	cerr << "Atomic network substitution aborted because O atom bonded to " 
	     << bonds[i].size() << " other atoms instead of 2" << "\n"
	     << "Occurred for atom " << i << "\n";
	return false;
      }	
    }
    else{
      cerr << "Atomic network substitution aborted because atom type other than "
              "Si or O detected" << "\n" << "Occurred for atom " << i << "\n";
      return false;
    }
  }

  int firstSiID = -1;
  for(int i = 0; i < numAtoms; i++){
    if(origNet->atoms[i].type.compare("Si") == 0){
      firstSiID = i;
      break;
    }
  }
  
  if(firstSiID == -1){
    cerr << "Error: Atom substitution failed because structure does not "
            "contain any Si atoms \n";
    return false;
  }

  vector<bool> atomProc = vector<bool>(numAtoms, false); // Whether atom has been processed
  vector<bool> atomSubs = vector<bool>(numAtoms, false); // Whether atom has been substituted
  int numProc = 0; // Number of processed atoms
     
  vector< pair<int, bool> > atomsToProcess; // Stack of atom/substituion info
  vector<int> fromIDs; // Stack of where each atom/subst came from

  pair<int, bool> seed(firstSiID, substituteSeed);
  atomsToProcess.push_back(seed);  // Seed with the first Si atom and not substituting it
  fromIDs.push_back(-1); // Starting node

  while(atomsToProcess.size() != 0){
    pair<int, bool> atomInfo = atomsToProcess.back(); atomsToProcess.pop_back();
    int id = atomInfo.first; bool subst = atomInfo.second;
    int fromID = fromIDs.back(); fromIDs.pop_back();

    // Atom previously processed. Ensure that atom substitution agrees with previous results for Si atoms
    if(atomProc[id] ){
      if((origNet->atoms[id].type.compare("Si") == 0) && (atomSubs[id] != subst)){
	cerr << "Atomic network substitution aborted because every other Si atom "
	        "substitution criteria failed" << "\n";
	return false;
      }
    }
    // First time atom is processed
    else {
      numProc++;
      // Mark atom as processed
      atomProc[id] = true;
      
      // If oxygen, should not be substituted but reverse substituion state for next set of Si atoms
      if(origNet->atoms[id].type.compare("O") == 0){
	atomSubs[id] = false;
	subst = !subst;
      }
      // If silicon, pass on substitution state to oxygen atoms, who will then reverse it
      // Record substitution status
      else if(origNet->atoms[id].type.compare("Si") == 0){
	atomSubs[id] = subst;
      }
      // Invalid type
      else {
	cerr << "Atomic network substitution aborted because atom type other than Si or O detected" << "\n"
	     << "Occurred for atom " << id << "\n";
	return false;
      }
	
      // Added bonded atoms to stack with appropriate subsitution state
      for(unsigned int j = 0; j < bonds[id].size(); j++){
	if(bonds[id][j] != fromID){
	  atomsToProcess.push_back(pair<int, bool>(bonds[id][j], subst));
	  fromIDs.push_back(id);
	}
      }
    }
  }

  if(numProc != numAtoms){
    cerr << "Atom network substituion failed because not all atoms are interconnected" << "\n"
	 << "Visited " << numProc << " out of " << numAtoms << " atoms " << "\n";
    return false;
  }

  // Perform final consistency check
  for(int i = 0; i < numAtoms; i++){
    // Just check pair of atoms bonded to each oxygen
    if(origNet->atoms[i].type.compare("O") == 0){
      if(atomSubs[bonds[i][0]] && atomSubs[bonds[i][1]]){ 
	cerr << "Atom network substitution failed final consistency check." << "\n"
	     << "Oxygen atom #" << i << " is bound to atoms " <<  bonds[i][0] << " and " << bonds[i][1] << " which are of same type. " << "\n"
	     << "Algorithmic issue likely." << "\n"
	     << "Please contact the source code provider" << "\n";
	return false;
      }
    }
  }

  // Copy original information
  origNet->copy(newNet);
  
  // Change atoms that should be subsituted
  *numSubstitutions = 0;
  for(int i = 0; i < numAtoms; i++){
    if(atomSubs[i]){
      newNet->atoms[i].type   = "Al";
      newNet->atoms[i].radius = lookupRadius("Al", radial);
      (*numSubstitutions)++;
    }
  }

  // Successful substitution
  return true;
}


ptrdiff_t myrandom (ptrdiff_t i) { return rand();} // Random number generator function for fractional substitution
ptrdiff_t (*p_myrandom)(ptrdiff_t) = myrandom;     // Pointer to it

bool comparePairs(pair<int,int> p_one, pair<int,int> p_two){
  return p_one.first < p_two.first;
}



/* Attempt to substitute the specified fraction of Si atom with an Al atom. ATOM_NETWORK may only consist of
*  Si and O atoms, where each Si atom must be bonded to exactly 4 oxygen atoms and each oxygen atom must
*  be bonded to exactly 2 Si atoms. Returns true iff the substitution was successful and stores the number of substitutions using the
*  provided reference. The provided boolean specifies whether the seeded silicon atom is substituted or not. The function works
*  by first substituting every other Si atom and then reverting some of the substituted atoms back to Si. The provided
*  random number generator seed is used to choose which atoms to switch back.*/
bool fracSubstituteAtoms(ATOM_NETWORK *origNet, ATOM_NETWORK *newNet, 
                         bool substituteSeed, double frac, int randSeed, 
                         int *numSubstitutions, double *fracSub, bool radial){
  int numAtoms = origNet->numAtoms;
  double max_bond_length = 1.95;

  if(frac > 0.5){
    cerr << "Fractional atomic network substitution failed because the fraction can not exceed 0.5" << "\n";
    return false;
  }
  if(frac < 0.0){
    cerr << "Fractional atomic network substitution failed because of invalid negative fraction " << frac << "\n";
    return false;
  }
  
  vector< vector<int> > bonds = vector< vector<int> >(numAtoms, vector<int>());
  for(int i = 0; i < numAtoms; i++){
    ATOM atom_one = origNet->atoms[i];
    for(int j = i + 1; j < numAtoms; j++){
      ATOM atom_two = origNet->atoms[j];
      if(origNet->calcDistance(atom_one, atom_two) < max_bond_length){
	if(atom_one.type.compare(atom_two.type) == 0){
	  cerr << "Fractional atomic network substitution aborted because atoms of same type are bonded to one another" << "\n"
	       << "Occurred for type " << atom_one.type << " between atoms " << i << " and " << j << "\n";
	  return false;
	}
	bonds[i].push_back(j);
	bonds[j].push_back(i);
      }									    
    }
  }

  int numSi = 0;

  for(int i = 0; i < numAtoms; i++){
    if(origNet->atoms[i].type.compare("Si") == 0){
      if(bonds[i].size() != 4){
	cerr << "Fractional atomic network substitution aborted because Si atom "
	        "bonded to " << bonds[i].size() << " other atoms instead of 4" << "\n"
	     << "Occurred for atom " << i << "\n";
	return false;
      }
      numSi++;
    }
    else if(origNet->atoms[i].type.compare("O") == 0){
      if(bonds[i].size() != 2){
	cerr << "Fractional atomic network substitution aborted because O atom "
	        "bonded to " << bonds[i].size() << " other atoms instead of 2" << "\n"
	     << "Occurred for atom " << i << "\n";
	return false;
      }	
    }
    else{
      cerr << "Fractional atomic network substitution aborted because atom type "
              "other than Si or O detected" << "\n"
	       << "Occurred for atom " << i << "\n";
      return false;
    }
  }

  int firstSiID = -1;
  for(int i = 0; i < numAtoms; i++){
    if(origNet->atoms[i].type.compare("Si") == 0){
      firstSiID = i;
      break;
    }
  }
  
  if(firstSiID == -1){
    cerr << "Error: Fractional atomic network substitution failed because "
            "structure does not contain any Si atoms \n";
    return false;
  }

  vector<bool> atomProc = vector<bool>(numAtoms, false); // Whether atom has been processed
  vector<bool> atomSubs = vector<bool>(numAtoms, false); // Whether atom has been substituted
  int numProc = 0; // Number of processed atoms
     
  vector< pair<int, bool> > atomsToProcess; // Stack of atom/substituion info
  vector<int> fromIDs; // Stack of where each atom/subst came from

  pair<int, bool> seed(firstSiID, substituteSeed);
  atomsToProcess.push_back(seed);  // Seed with the first Si atom and not substituting it
  fromIDs.push_back(-1); // Starting node

  while(atomsToProcess.size() != 0){
    pair<int, bool> atomInfo = atomsToProcess.back(); atomsToProcess.pop_back();
    int id = atomInfo.first; bool subst = atomInfo.second;
    int fromID = fromIDs.back(); fromIDs.pop_back();

    // Atom previously processed. 
    // Ensure that atom substitution agrees with previous results for Si atoms
    if(atomProc[id] ){
      if((origNet->atoms[id].type.compare("Si") == 0) && (atomSubs[id] != subst)){
	cerr << "Fractional atomic network substitution aborted because every other "
	        "Si atom substitution criteria failed" << "\n";
	return false;
      }
    }
    // First time atom is processed
    else {
      numProc++;
      // Mark atom as processed
      atomProc[id] = true;
      
      // If oxygen, should not be substituted but reverse substituion state for next set of Si atoms
      if(origNet->atoms[id].type.compare("O") == 0){
	atomSubs[id] = false;
	subst = !subst;
      }
      // If silicon, pass on substitution state to oxygen atoms, who will then reverse it
      // Record substitution status
      else if(origNet->atoms[id].type.compare("Si") == 0){
	atomSubs[id] = subst;
      }
      // Invalid type
      else {
	cerr << "Fractional atomic network substitution aborted because atom type "
	        "other than Si or O detected" << "\n"
	     << "Occurred for atom " << id << "\n";
	return false;
      }
	
      // Added bonded atoms to stack with appropriate subsitution state
      for(unsigned int j = 0; j < bonds[id].size(); j++){
	if(bonds[id][j] != fromID){
	  atomsToProcess.push_back(pair<int, bool>(bonds[id][j], subst));
	  fromIDs.push_back(id);
	}
      }
    }
  }

  if(numProc != numAtoms){
    cerr << "Fractional atom network substituion failed because not all atoms "
            "are interconnected" << "\n" << "Visited " << numProc << " out of " 
         << numAtoms << " atoms " << "\n";
    return false;
  }

  // Perform final consistency check
  for(int i = 0; i < numAtoms; i++){
    // Just check pair of atoms bonded to each oxygen
    if(origNet->atoms[i].type.compare("O") == 0){
      if(atomSubs[bonds[i][0]] && atomSubs[bonds[i][1]]){ 
	cerr << "Fractional atom network substitution failed final consistency check." 
	     << "\n" << "Oxygen atom #" << i << " is bound to atoms " <<  bonds[i][0] 
	     << " and " << bonds[i][1] << " which are of same type. " << "\n"
	     << "Algorithmic issue likely." << "\n"
	     << "Please contact the source code provider" << "\n";
	return false;
      }
    }
  }
  
  // Copy original information
  origNet->copy(newNet);

  int totalSubs = nearestInt(frac*numSi); // How many substituions to actually perform
  *fracSub = (1.0*totalSubs)/numSi;        // Store the actual substitution fraction
  
  srand(randSeed); // Seed random number generator

  // Record all originally substituted atoms
  vector<pair<int, int> > subIDs;
  for(int i = 0; i < numAtoms; i++){
    if(atomSubs[i]){
      subIDs.push_back(pair<int,int> (rand(), i));
    }
  }

  sort(subIDs.begin(), subIDs.end(), comparePairs); // Sort the pairs by their random values
 
  // Change atoms that should be substituted
  for(*numSubstitutions = 0; *numSubstitutions < totalSubs; (*numSubstitutions)++){
    int id = subIDs[*numSubstitutions].second;
    newNet->atoms[id].type   = "Al";
    newNet->atoms[id].radius = lookupRadius("Al", radial);
  }


  // Change Oxygen type if O atom connected to Al
  for(int i=0; i < numAtoms; i++){
    if(origNet->atoms[i].type.compare("O") == 0)
      {
      if(newNet->atoms[(bonds[i].at(0))].type.compare("Al") == 0 || newNet->atoms[(
         bonds[i].at(1))].type.compare("Al") == 0)
        {
        newNet->atoms[i].type   = "O_Al";
        };
      };  
    };


 
  // Successful substitution
  return true;
}


/* Maciek's substitution function (based on Thomas -fsub function) */
/* Attempt to substitute the specified fraction of Si atom with an Al atom. 
 * ATOM_NETWORK may only consist of Si and O atoms, where each Si atom must be 
 * bonded to exactly 4 oxygen atoms and each oxygen atom must be bonded to 
 * exactly 2 Si atoms. Returns true if the substitution was successful and 
 * stores the number of substitutions using the provided reference. The 
 * provided boolean specifies whether the seeded silicon atom is substituted
 * or not. This function does not require initial 50/50 substitution as in 
 * Thomas version, and therefore is suitable for systems in which not every 
 * other atom is a Si atom */
bool fracSubstituteAtoms_Maciek(ATOM_NETWORK &origNet, ATOM_NETWORK &newNet, 
                                bool substituteSeed, double frac, int randSeed, 
                                int &numSubstitutions, double &fracSub, bool radial){
  int numAtoms = origNet.numAtoms;
  double max_bond_length = 1.95;

  if(frac > 0.50){
    cerr << "Fractional atomic network substitution failed because the fraction "
            "can not exceed 0.5" << "\n";
    return false;
  }
  if(frac < 0.0){
    cerr << "Fractional atomic network substitution failed because of invalid "
            "negative fraction " << frac << "\n";
    return false;
  }
  
  vector< vector<int> > bonds = vector< vector<int> >(numAtoms, vector<int>());
  for(int i = 0; i < numAtoms; i++){
    ATOM atom_one = origNet.atoms[i];
    for(int j = i + 1; j < numAtoms; j++){
      ATOM atom_two = origNet.atoms[j];
      if(origNet.calcDistance(atom_one, atom_two) < max_bond_length){
	if(atom_one.type.compare(atom_two.type) == 0){
	  cerr << "Fractional atomic network substitution aborted because atoms of "
	          "same type are bonded to one another" << "\n" << "Occurred for type "
	       << atom_one.type << " between atoms " << i << " and " << j << "\n";
	  return false;
	}
	bonds[i].push_back(j);
	bonds[j].push_back(i);
      }									    
    }
  }

  int numSi = 0;

  for(int i = 0; i < numAtoms; i++){
    if(origNet.atoms[i].type.compare("Si") == 0){
      if(bonds[i].size() != 4){
	cerr << "Fractional atomic network substitution aborted because Si atom "
	        "bonded to " << bonds[i].size() << " other atoms instead of 4" << "\n"
	     << "Occurred for atom " << i << "\n";
	return false;
      }
      numSi++;
    }
    else if(origNet.atoms[i].type.compare("O") == 0){
      if(bonds[i].size() != 2){
	cerr << "Fractional atomic network substitution aborted because O atom "
	        "bonded to " << bonds[i].size() << " other atoms instead of 2" << "\n"
	     << "Occurred for atom " << i << "\n";
	return false;
      }	
    }
    else{
      cerr << "Fractional atomic network substitution aborted because atom type "
              "other than Si or O detected" << "\n"
	       << "Occurred for atom " << i << "\n";
      return false;
    }
  } // analysis of the input structure completed

  // save IDs of Si atoms in SiatomsIDs vector
  vector <int> SiatomsIDs;
  for(int i = 0; i < numAtoms; i++){
    if(origNet.atoms[i].type.compare("Si") == 0){
      SiatomsIDs.push_back(i);
    }
  }


  int totalSubs = nearestInt(frac*numSi); // How many substituions to actually perform
  fracSub = (1.0*totalSubs)/numSi;        // Store the actual substitution fraction
  
  srand(randSeed); // Seed random number generator

  // Record all Si atoms (each of them can be substituted) 
  vector<pair<int, int> > subIDs;
  for(int i = 0; i < numSi; i++){
      subIDs.push_back(pair<int,int> (rand(), i));
    }
  
  sort(subIDs.begin(), subIDs.end(), comparePairs); // Sort the pairs by their random values

  vector<bool> atomSubs = vector<bool>(numAtoms, false); // Whether atom has been substituted
 
  // Change atoms that should be substituted
  for(numSubstitutions = 0; numSubstitutions < totalSubs; numSubstitutions++){
    atomSubs[SiatomsIDs[subIDs[numSubstitutions].second]]=true;
  } // now atomsSubs points all Si atoms substituted in the intial step 


  // Verify if the random distribution of Al atoms satisfy the condition not to have two neighboring Al atoms (Lowenstein's rule)
  // Fix if neccessary
  int counter=0;
  int problematicAl=-1;

  int max_try=10000; // number of swaps when trying to find a configuration satisfying Lowenstein's rule

  do{
  problematicAl=-1;
  // Step 1: identify atom with problems    
  for(int i = 0; i < numAtoms; i++){
    if(origNet.atoms[i].type.compare("Si") == 0 &&atomSubs[i]==true)
      {
      // find list of 4 neighboring Si/Al atoms
      vector <int> nlist;
      for(int j=0;j<4;j++)
         {
         if(bonds[bonds[i][j]][0]==i) {nlist.push_back(bonds[bonds[i][j]][1]);} else {nlist.push_back(bonds[bonds[i][j]][0]);};
         }
      
      if(atomSubs[nlist[0]]||atomSubs[nlist[1]]||atomSubs[nlist[2]]||atomSubs[nlist[3]]) // if any of neighbors is Al, move Al from position i
        {
        problematicAl=i; // Al atom i needs to be this shifted;
        break;
        }
      }
  }

  // Step 2: identify a suitable place for Al atom
  int suitableSi=-1;
  if(problematicAl>-1)
  {

    for(int i = 0; i < numAtoms; i++){
      if(origNet.atoms[i].type.compare("Si") == 0 && atomSubs[i]==false)
      {
        // find list of 4 neighboring Si/Al atoms
        vector <int> nlist;
        for(int j=0;j<4;j++)
        {
           if(bonds[bonds[i][j]][0]==i) {
	         nlist.push_back(bonds[bonds[i][j]][1]);
	       } else {
		     nlist.push_back(bonds[bonds[i][j]][0]);
		   };
        }
      
        if(!atomSubs[nlist[0]] && !atomSubs[nlist[1]] && 
           !atomSubs[nlist[2]]&&!atomSubs[nlist[3]]) // if all of neighbors are Si, position i can be used
        {
          suitableSi=i; // Si atom i can be replaced with Al;
          break;
        }
      }
    }
  } // end of Step 2
    
  // Step 3: Swap 1 & 2
  if(problematicAl>-1&&suitableSi>-1)
    {
    atomSubs[problematicAl]=false;
    atomSubs[suitableSi]=true;
    }else{

     if(problematicAl>-1) 
       {
       cerr << "Cannot fix a random Al distribution by swapping Al atoms" << "\n";
       return false;
       }
    }
  // Step 4: increase counter
  counter++;
  }while(problematicAl>-1&&counter<max_try*numSi);


  if(counter==(max_try*numSi))
    {
    cerr << "Could not fix the initial Al distribution in " << max_try*numSi << "steps\n";
    return false;
    }



  // Now update the original structure according to the generate distribution of Al atoms
  
  // Copy original information
  origNet.copy(&newNet);


  // Change atoms that should be substituted
  for(int i = 0; i < numAtoms; i++){
    if(atomSubs[i]){
      newNet.atoms[i].type   = "Al";
      newNet.atoms[i].radius = lookupRadius("Al", radial);
    }
  }


 
  // Change Oxygen type if O atom connected to Al
  for(int i=0; i < numAtoms; i++){
    if(origNet.atoms[i].type.compare("O") == 0)
      {
      if(newNet.atoms[(bonds[i].at(0))].type.compare("Al") == 0 || newNet.atoms[(bonds[i].at(1))].type.compare("Al") == 0)
        {
        newNet.atoms[i].type   = "O_Al";
        };
      };  
    };

 
  // Successful substitution
  return true;
}





/* Returns the integer nearest to the provided double.*/
int nearestInt(double num){
  return (int)(floor(num + 0.5));
}

/* Determines whether a specific supercell size satisfies the non-overlapping sphere requirement */
int check_sphere_overlap(int num_a, int num_b, int num_c, double diam, ATOM_NETWORK *atmnet) {
  //check each neighbouring cell in three dimensions and find shortest distance to image
  //it is sufficient to check only those images (a b c) where the leftmost non-zero term is one, e.g. (0 1 -1) but we can skip (0 -1 1), since it is equivalent in magnitude
  double min_d = -1.0; //set min distance to invalid quantity - if we cannot find any result we can return the error case
  bool overlaps = false;
  for(int a=0; a<=1 && !overlaps; a++) {
    for(int b=-1; b<=1 && !overlaps; b++) {
      for(int c=-1; c<=1 && !overlaps; c++) {
        if (
        (a==0 && b==0 && c==1)
        ||
        (a==0 && b==1)
        ||
        (a==1)
        )
        {
          Point image_abc(num_a*a, num_b*b, num_c*c);
          Point image_xyz = atmnet->abc_to_xyz(image_abc);
          double d=image_xyz.magnitude();
          if(d<min_d || min_d<0) {
            min_d = d;
            if(min_d<diam+0.001) overlaps = true;
          }
        }
      }
    }
  }
  if(min_d<0) return -1; //bad_angles!
  else if(overlaps) return 0; //wrong
  else return 1; //correct
}

/* Returns the id of the VOR_NODE in the provided VORONOI_NETWORK whose coordinates
 * match those of the provided Point.
 */
int getNodeID(Point pt, ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet){
  double minDist = DBL_MAX;
  int minID = -1;
  for(unsigned int i = 0; i < vornet->nodes.size(); i++){
    VOR_NODE curNode = vornet->nodes.at(i);
    double dist = atmnet->calcDistanceXYZ(pt[0], pt[1], pt[2], curNode.x, curNode.y, curNode.z);
    if(dist < threshold)
      return i;
    else {
      if(dist < minDist){
	minDist = dist;
	minID = i;
      }
    }
  }
 
  cerr << "Warning : When identifying Voronoi node, the distance exceeded the threshold of " << threshold << "\n"
       << "Occurred during analysis of " << atmnet->name << "\n"
       << "Closest node was within " << minDist << "\n"
       << "Proceeding with analysis" << "\n";
  return minID;
}

