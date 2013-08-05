/** This file contains all of the necessary geometric object classes and
 *  procedures that are used throughout the package.
 */

/** Simple class used to store a set of three numbers and perform 
 *  some trivial geometric operations. */

#include <cstdlib>
#include <cmath>
#include <float.h>      //change DBL_MAX to std::numeric_limits<double>::max in future from #include <limits>

#include "zeo_consts.h"
#include "geometry.h"

using namespace std;

Point::Point (double myX, double myY, double myZ){  
  vals[0] = myX; 
  vals[1] = myY; 
  vals[2] = myZ;
}
const Point Point::add(Point other) const { 
  return Point(vals[0] + other[0], vals[1] + other[1], vals[2] + other[2]); 
}
const bool Point::equals(Point other) const { 
  double tol = 0.00000001; 
  return abs(vals[0] - other[0]) < tol && abs(vals[1] - other[1]) < tol && abs(vals[2] - other[2]) < tol; 
}
const Point  Point::subtract(Point other) const { 
  return Point(vals[0]-other[0], vals[1]-other[1], vals[2]-other[2]); 
}
const double Point::magnitude() const {
  return sqrt(vals[0]*vals[0] + vals[1]*vals[1] + vals[2]*vals[2]); 
}
const Point  Point::scale(const double& factor) const { 
  return Point(vals[0]*factor, vals[1]*factor, vals[2]*factor); 
}
const double Point::dot_product(Point other) const {
  return vals[0]*other[0] + vals[1]*other[1] + vals[2]*other[2]; 
}
const Point Point::unit() const {
  return Point(vals[0]/magnitude(),vals[1]/magnitude(),vals[2]/magnitude());
}
const Point Point::cross(Point other) const {
  return Point(vals[1]*other[2]-vals[2]*other[1],vals[2]*other[0]-vals[0]*other[2],vals[0]*other[1]-vals[1]*other[0]);
}
double& Point::operator[](const int& index) {
     if(index < 0 || index > 2){
      cerr << "Error: Invalid index to [] operator for Point instance" << "\n"
	   << "Exiting..." << "\n";
      exit(1);
    }
    return vals[index];
}
const Point Point::operator-(Point other) const {
   return Point(vals[0]-other[0], vals[1]-other[1], vals[2]-other[2]); 
}
const Point Point::operator+(Point other) const {
   return Point(vals[0] + other[0], vals[1] + other[1], vals[2] + other[2]); 
}
const double Point::operator*(Point other) const {
  return vals[0]*other[0] + vals[1]*other[1] + vals[2]*other[2]; 
}
void Point::print(std::ostream &out) const {
  out << vals[0] << " " << vals[1] << " " << vals[2];
}
ostream &operator<<(ostream &out, Point &obj) {
    obj.print(out);
}


bool pointIsLess(Point p1, Point p2){
  double tol = 0.0000001;
  if (abs(p1[0] - p2[0]) > tol)
    return p1[0] < p2[0];
  else if (abs(p1[1] - p2[1]) > tol)
    return p1[1] < p2[1];
  else if (abs(p1[2] - p2[2]) > tol)
    return p1[2] < p2[2];
  else
    return false;
}


/* Construct a TRIPLET with the three provided components.*/
TRIPLET::TRIPLET(int myX, int myY, int myZ){
  vals[0] = x = myX; 
  vals[1] = y = myY; 
  vals[2] = z = myZ;
}

/* Access one of the TRIPLETs three values.*/
int& TRIPLET::operator[](int index){
  if(index < 0 || index > 2){
    cerr << "Error: Invalid index to [] operator for TRIPLET instance" << "\n"
	 << "Exiting..." << "\n";
    exit(1);
  }
  return vals[index];
}

/** Add each component of the triplet to that of the provided TRIPLET
 *  and return the result. */
const TRIPLET TRIPLET::operator+(const TRIPLET& other){
  return TRIPLET(x + other.x, y + other.y, z + other.z);
}


XYZ::XYZ(double myX, double myY, double myZ) : x(myX), y(myY), z(myZ) {}
void XYZ::print(ostream &out) const {
  cout << "x:" << x << "   y:" << y << "   z:" << z << "\n";
}
const double XYZ::magnitude() const {
  return sqrt(x*x+y*y+z*z);
}

/** Return a instance whose three components are equal to this XYZ's components
 *  multiplied by the provided factor. */
const XYZ XYZ::scale(const double& factor) const{
  return XYZ(factor*x, factor*y, factor*z);
}
const XYZ XYZ::unit() const {
  return XYZ(x/magnitude(),y/magnitude(),z/magnitude());
}
const XYZ XYZ::cross(const XYZ& other) const{
  return XYZ(y*other.z-z*other.y,z*other.x-x*other.z,x*other.y-y*other.x);
}
const double XYZ::dot_product(const XYZ& other) const {
  return x*other.x + y*other.y + z*other.z; 
}
double XYZ::angle_between(const XYZ& other)  const {
  double cos_angle = dot_product(other)/(magnitude()*other.magnitude());
  // These steps are necessary for rounding issues when cos_angle should be 
  // exactly 1 or -1, but falls very slightly on the wrong side
  if(cos_angle>1) {
    cos_angle = 1;
  } else if(cos_angle<-1) {
    cos_angle = -1; 
  }
  double output = acos(cos_angle);
  if(isnan(output)) {
    return 0;
  } else return output; //this step is necessary in the case of identity, which produces NAN
}
double& XYZ::operator[](const int& index){
     if (index==0)
       return x;
     if (index==1)
       return y;
     if (index==2)
       return z;
     cerr << "Error: Invalid index to [] operator for XYZ instance" << "\n"
	  << "Exiting..." << "\n";
     exit(1);
}
const XYZ XYZ::operator-(const XYZ& other) const {
   return XYZ(x-other.x, y-other.y, z-other.z); 
}
const XYZ XYZ::operator+(const XYZ& other) const {
   return XYZ(x + other.x, y + other.y, z + other.z); 
}
const double XYZ::operator*(const XYZ& other) const {  //Dot product
  return x*other.x + y*other.y + z*other.z; 
}




/** Convert coordinates relative to the unit cell (a,b,c) to absolute cartesian
    coordinates (x,y,z). */
Point abc_to_xyz(double a, double b, double c, const XYZ& v_a, const XYZ& v_b, const XYZ& v_c){

  double x = v_a.x+b*v_b.x+c*v_c.x;
  double y = b*v_b.y+c*v_c.y;
  double z = c*v_c.z;
  return Point (x, y, z);
}



double calcSphereVolume(double radius){
  return 4.0/3.0*PI*radius*radius*radius;
}

/** Calculates the Euclidean distance between (x1,y1,z1) and (x2,y2,z2). */
double calcEuclideanDistance(double x1, double y1, double z1, double x2, double y2, double z2){
  double dx = x1 - x2;  double dy = y1 - y2; double dz = z1 - z2;
  return sqrt(dx*dx + dy*dy + dz*dz);
}

/** Calculates the Euclidean distance between the provided Points. */
double calcEuclideanDistance(Point p1, Point p2){
  return calcEuclideanDistance(p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]);
}







/** Calculates the number of intersections between the line that travels through Points p1 and p2
 *  and the sphere with the provided characterstics. Returns a pair made up of the number of intersections
 *  and a vector containing the intersection points. Refer to http://local.wasp.uwa.edu.au/~pbourke/geometry/sphereline
 *  for more information.
 */
pair<int, vector<Point> > findLineIntersections(Point p1, Point p2, Point circCenter, double rad){

  double tolerance = 0.00001; // Tolerance for determining whether two numbers are equal

  double a = (p2[0] - p1[0])*(p2[0] - p1[0]) + (p2[1] - p1[1])*(p2[1] - p1[1]) + (p2[2] - p1[2])*(p2[2] - p1[2]);
  double b = 2.0*((p2[0] - p1[0])*(p1[0] - circCenter[0]) + (p2[1] - p1[1])*(p1[1] - circCenter[1]) + (p2[2] - p1[2])*(p1[2] - circCenter[2]));
  double c = circCenter[0]*circCenter[0] + circCenter[1]*circCenter[1] + circCenter[2]*circCenter[2] 
    + p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2] 
    - 2.0*(circCenter[0]*p1[0] +circCenter[1]*p1[1] + circCenter[2]*p1[2]) - rad*rad;

  double discriminant = b*b - 4.0*a*c;

  int numResults;
  vector<Point> results;

  if(abs(discriminant) < tolerance){
    // Discriminant is zero. Line is tangent and therefore only has one intersection
    numResults = 1;
    double u = (-b/(2.0*a));
    Point intersection = p1.add((p2.subtract(p1)).scale(u));
    results.push_back(intersection);
  }
  else if (discriminant < 0){
    //Disciminant is negative. No intersection between line and sphere.
    numResults = 0;
  }
  else{
    //Discriminant is positive. Two intersection points exist
    numResults = 2;

    double u1 = ((-b + sqrt(discriminant))/(2.0*a));
    Point intersection1 = p1.add((p2.subtract(p1)).scale(u1));
    results.push_back(intersection1);

    double u2 = ((-b - sqrt(discriminant))/(2.0*a));
    Point intersection2 = p1.add((p2.subtract(p1)).scale(u2));
    results.push_back(intersection2);
  }
  return pair<int, vector<Point> > (numResults, results);
}


pair<int, vector<Point> > findLineSegIntersections(Point p1, Point p2, Point circCenter, double rad){
  pair<int, vector<Point> > prelimResults = findLineIntersections(p1, p2, circCenter, rad);
 
  if(prelimResults.first == 0)
    return prelimResults;
  else{
    double tolerance = 0.00001;
    int numResults = 0;
    vector<Point> pts = vector<Point> ();
    double lineSegLength = calcEuclideanDistance(p1, p2);
    for(int i = 0; i < prelimResults.first; i++){
      Point possSol = prelimResults.second.at(i);
      double p1Dist = calcEuclideanDistance(possSol, p1);
      double p2Dist = calcEuclideanDistance(possSol, p2);
      if((p1Dist <= lineSegLength + tolerance) && (p2Dist <= lineSegLength + tolerance)){
	pts.push_back(possSol);
	numResults++;
      }
    }
    return pair<int, vector<Point> > (numResults, pts);
  }
}

/* Returns the shortest distance along a sphere of the provided circle radius 
 * between the two provided points. Refer to http://en.wikipedia.org/wiki/Great-circle_distance. */ 
double calcSphereDistance(Point p1, Point p2, double circRad){
  pair<double,double> coords1 = findLongAndLat(p1);
  pair<double,double> coords2 = findLongAndLat(p2);
  double phi1 = coords1.first; double phi2 = coords2.first;
  double lambda1 = coords1.second; double lambda2 = coords2.second;
  double deltaLambda = lambda1-lambda2;

  double deltaSigma = atan( (sqrt (pow(cos(phi1)*sin(deltaLambda),2.0) +  pow(cos(phi2)*sin(phi1)-sin(phi2)*cos(phi1)*cos(deltaLambda),2))) 
			    / (sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(deltaLambda)));
  return circRad*deltaSigma;
}

/* Returns the projection of the provided Point onto the circle with the specified center and radius. */
Point projectPointOnSphere(Point p, double circRad, Point circCenter){
  double deltaX = p[0] - circCenter[0];
  double deltaY = p[1] - circCenter[1];
  double deltaZ = p[2] - circCenter[2];
  double factor = sqrt((circRad*circRad)/(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ));
  return Point(factor*deltaX+circCenter[0], factor*deltaY+circCenter[1], factor*deltaZ+circCenter[2]);
}

/* Returns a pair of doubles representing the longitude and latitude of the given (x,y,z) point.
 * Refer to http://en.wikipedia.org/wiki/Geodetic_system#geodetic_to.2Ffrom_ECEF_coordinates */
pair<double,double> findLongAndLat(Point pt){
  double r = sqrt(pt[0]*pt[0] + pt[1]*pt[1]);
  double phi = atan(pt[2]/r);
  double lambda = atan2(pt[1],pt[0]);
  return pair<double,double> (lambda, phi);
}

/** Returns the shortest distance from point to plane **/
double distToPlane(Point pnt,Point p,Point normal){
  if ((normal * (pnt - p)) < 0)
    return -(normal * (pnt - p));
  
  return normal * (pnt - p);
}

