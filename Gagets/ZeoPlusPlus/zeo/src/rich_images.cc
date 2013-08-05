#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <queue>
//Rich edit
//#include "cell.cc"

typedef struct {int x; int y; int z;} coords_int; //standard triplet of ints
//end Rich edit

const int max_shells(10);
const int ms2(2*max_shells+1);

double bx,bxy,by,bxz,byz,bz;
double bxinv,byinv,bzinv,ivol;

voronoicell unitcell,c;

bool unit_cell_test(int i,int j,int k) {
	double x=i*bx+j*bxy+k*bxz,y=j*by+k*byz,z=k*bz;
	double rsq=x*x+y*y+z*z;
	return unitcell.plane_intersects(x,y,z,rsq);
}

bool unit_cell_intersect(int l) {
	int i,j;
	if(unit_cell_test(l,0,0)) return true;
	for(i=1;i<l;i++) {
		if(unit_cell_test(l,i,0)) return true;
		if(unit_cell_test(-l,i,0)) return true;
	}
	for(i=-l;i<=l;i++) if(unit_cell_test(i,l,0)) return true;
	for(i=1;i<l;i++) for(j=-l+1;j<=l;j++) {
		if(unit_cell_test(l,j,i)) return true;
		if(unit_cell_test(-j,l,i)) return true;
		if(unit_cell_test(-l,-j,i)) return true;
		if(unit_cell_test(j,-l,i)) return true;
	}
	for(i=-l;i<=l;i++) for(j=-l;j<=l;j++) if(unit_cell_test(i,j,l)) return true;
	return false;
}

// To carry out the computation centered on the origin, rather than the domain
// center, change the 1's in the plane calls to 2,0,2,0,2,0.
bool intersects(int i,int j,int k,double &vol) {
	c.init(unitcell);
	c.translate(i*bx+j*bxy+k*bxz,j*by+k*byz,k*bz);
	if(!c.plane(0,0,bzinv,1)) return false;
	if(!c.plane(0,0,-bzinv,1)) return false;
	if(!c.plane(0,byinv,-byz*byinv*bzinv,1)) return false;
	if(!c.plane(0,-byinv,byz*byinv*bzinv,1)) return false;
	if(!c.plane(bxinv,-bxy*bxinv*byinv,(bxy*byz-by*bxz)*ivol,1)) return false;
	if(!c.plane(-bxinv,bxy*bxinv*byinv,(-bxy*byz+by*bxz)*ivol,1)) return false;
	vol=c.volume()*ivol;
	return true;
}

inline void unit_cell_apply(int i,int j,int k) {
	double x=i*bx+j*bxy+k*bxz,y=j*by+k*byz,z=k*bz;
	unitcell.plane(x,y,z);
	unitcell.plane(-x,-y,-z);
}

void compute_unit_cell() {
	int i,j,l=1;

	// Initialize the Voronoi cell to be a very large rectangular box
	const double ucx=max_shells*bx,ucy=max_shells*by,ucz=max_shells*bz;
	unitcell.init(-ucx,ucx,-ucy,ucy,-ucz,ucz);

	// Repeatedly cut the cell by shells of periodic image particles
	while(l<2*max_uc_shells) {

		// Check to see if any of the planes from the current shell
		// will cut the cell
		if(unit_cell_intersect(l)) {

			// If they do, apply the plane cuts from the current
			// shell
			unit_cell_apply(l,0,0);
			for(i=1;i<l;i++) {
				unit_cell_apply(l,i,0);
				unit_cell_apply(-l,i,0);
			}
			for(i=-l;i<=l;i++) unit_cell_apply(i,l,0);
			for(i=1;i<l;i++) for(j=-l+1;j<=l;j++) {
				unit_cell_apply(l,j,i);
				unit_cell_apply(-j,l,i);
				unit_cell_apply(-l,-j,i);
				unit_cell_apply(j,-l,i);
			}
			for(i=-l;i<=l;i++) for(j=-l;j<=l;j++) unit_cell_apply(i,j,l);
		} else return;
		l++;
	}

	// If the routine makes it here, then the unit cell still hasn't been
	// completely bounded by the plane cuts. Give the memory error code,
	// because this is mainly a case of hitting a safe limit, than any
	// inherent problem.
	voropp_fatal_error("Periodic cell computation failed",VOROPP_MEMORY_ERROR);
}

//Rich edit: changing main() to do_images() which can be called by my main
vector<coords_int> do_images(int bx_, int bxy_, int by_, int bxz_, int byz_, int bz_) {
//end Rich edit

	int i,j,k;
	double vol;
//Rich edit: handling arguments which aren't in argv format etc
/*
	// Check the command line syntax
	if(argc!=7) {
		fprintf(stderr,"Syntax: ./images bx bxy by bxz byz bz\n");
		return VOROPP_CMD_LINE_ERROR;
	}
*/
	// Decode the command line arguments
/*
	bx=atof(argv[1]);
	bxy=atof(argv[2]);
	by=atof(argv[3]);
	bxz=atof(argv[4]);
	byz=atof(argv[5]);
	bz=atof(argv[6]);
*/
	bx=bx_;
	bxy=bxy_;
	by=by_;
	bxz=bxz_;
	byz=byz_;
	bz=bz_;
//end Rich edit

	// Create inverse lengths and inverse volume
	bxinv=1/bx;byinv=1/by;bzinv=1/bz;ivol=bxinv*byinv*bzinv;

	compute_unit_cell();
	bool a[ms2*ms2*ms2],*ac(a+max_shells*(1+ms2*(1+ms2))),*acp;
	for(i=0;i<ms2*ms2*ms2;i++) a[i]=true;

	*ac=false;
	queue<int> q;
	q.push(0);q.push(0);q.push(0);

//Rich edit: keeping track of cells
	vector<coords_int> images;
//end Rich edit
	while(!q.empty()) {
		i=q.front();q.pop();
		j=q.front();q.pop();
		k=q.front();q.pop();
		if(intersects(i,j,k,vol)) {
//Rich edit: counter
//			printf("%d %d %d %g\n",i,j,k,vol);
			coords_int temp = {i, j, k};
			images.push_back(temp);
//end Rich edit
			acp=ac+i+ms2*(j+ms2*k);
			if(k>-max_shells&&*(acp-ms2*ms2)) {q.push(i);q.push(j);q.push(k-1);*(acp-ms2*ms2)=false;}
			if(j>-max_shells&&*(acp-ms2)) {q.push(i);q.push(j-1);q.push(k);*(acp-ms2)=false;}
			if(i>-max_shells&&*(acp-1)) {q.push(i-1);q.push(j);q.push(k);*(acp-1)=false;}
			if(i<max_shells&&*(acp+1)) {q.push(i+1);q.push(j);q.push(k);*(acp+1)=false;}
			if(j<max_shells&&*(acp+ms2)) {q.push(i);q.push(j+1);q.push(k);*(acp+ms2)=false;}
			if(k<max_shells&&*(acp+ms2*ms2)) {q.push(i);q.push(j);q.push(k+1);*(acp+ms2*ms2)=false;}
		}
	}
//Rich edit: return vector
	return images;
//end Rich edit
}

