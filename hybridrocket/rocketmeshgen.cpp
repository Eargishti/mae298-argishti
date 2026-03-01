#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <functional>
#include "headerg.cpp"

typedef std::function<scalar(scalar, scalar)> radfunc;

#define SQRT_2 1.41421356237L

std::string Edgar = "hi";

struct Radial_Functions {

	radfunc r_inner;
	radfunc r_outer;
	double z_length = 0;
	void something(){
		long double *a = new long double[3];

//	MeshData<long double> MyMesh(a, "points");
	};
};


struct Circular_let {
	std::function<double(double phi)> r_inner;
	double z;
};


//MeshData<long double> M1();


class Boundary_Surface {
public:

	enum class Shape {
		CIRCULAR_PLATE,
		CYLINDRICAL_SURFACE,
	};

	Boundary_Surface(Radial_Functions p) : R1(p) {
		shape = Shape::CYLINDRICAL_SURFACE;
	}
	
	void set_function(radfunc Boundary_Surface::*target, radfunc source){
		this->*target = source;
	};


private:
	Radial_Functions R1;
	Circular_let C1;
	std::function<double(double zeta)>  phi_min;
	std::function<double(double zeta)>  phi_max;
	long double z_min, z_max;
	Shape shape;

	 
};


class Mesh{
public:
	Circular_let Inlet;
	Circular_let Outlet;
	Points points;
	Faces faces;
	Cells cells;
	Neighbour neighbour;
	Owner owner;

};

long double Calculate_rp_from_fineness(int mp, long double R){


return 0;
};


void ButterFlyMeshPoint(long double R, long double rp, int *idx, Points *points){

	
	long double R_square = 0.6 * R / SQRT_2;
//2 Iter rp - rp = R_square
//R_square = rp (2 Iter - 1)
	int Iter = (int) ( 0.5 + ( (0.5* R_square ) / rp) );
	for (int j = 0; j < Iter; j++){
		for (int i = 0; i < 2 * j + 1; i++){
			points->x[*idx] = rp * (long double) (2 * j + 1 - i);
			points->y[*idx] = rp * (long double) (i);
			(*idx)++;
		};

		for (int i = 0; i < 2 * j + 1; i++){
			points->x[*idx] = -rp * (long double) (i);
			points->y[*idx] = rp * (long double) (2 * j + 1 - i);
			(*idx)++;
		}
		for (int i = 0; i < 2 * j + 1; i++){
			points->x[*idx] = -rp * (long double) (2 * j + 1 - i);
			points->y[*idx] = -rp * (long double) (i);
			(*idx)++;
		}
		for (int i = 0; i < 2 * j + 1; i++){
			points->x[*idx] = rp * (long double) (i);
			points->y[*idx] = -rp * (long double) (2 * j + 1 - i);
			(*idx)++;
		}
	}

	long double dtheta = 2 * M_PIl / (4 * (2 * Iter -1) );
	long double theta = 0.0L;

	long double x0 = 0, y0 = 0;
	int idx2 = *idx;
	long double *xref = (long double *) malloc( 4 * (2 * Iter - 1));
	long double *yref = (long double *) malloc( 4 * (2 * Iter - 1));

    x0 = rp * (long double) (2 * Iter - 1) +rp * sqrtl(2.0L);
	y0 = 0.0L;
	points->x[*idx] = x0;
	points->y[*idx] = y0;
	xref[idx2] = x0;
	yref[idx2] = y0;
    
	(*idx)++;
	idx2++;
	theta += dtheta;

	for (int i = 1; i < 2 * Iter - 1; i++){

		long double thetap = theta - dtheta;

		long double m1 = (sinl(theta) -sinl(thetap)) / (cosl(theta) - cosl(thetap));
		long double m2 = (R * sinl(theta) - rp * (long double) (i) ) / (R * cosl(theta) - rp * (long double)(2 * Iter - 1 -i));

        long double b1 = y0 - m1 * x0;
		long double b2 = R * sinl(theta) - m2 * R * cosl(theta);

        points->x[*idx] = (b2 - b1) / (m1 - m2);
		points->y[*idx] = m1 * points->x[*idx] + y0 - m1* x0;
		xref[idx2] = points->x[*idx];
		yref[idx2] = points->y[*idx];
        x0 = points->x[*idx];
        y0 = points->y[*idx];
		(*idx)++;
		idx2++;
		theta += dtheta;
}

x0 = 0.0L;
y0 = rp * (long double) (2 * Iter - 1) + rp * sqrtl(2.0L);
	points->x[*idx] = x0;
	points->y[*idx] = y0;
	xref[idx2] = x0;
	yref[idx2] = y0;
    
	(*idx)++;
	idx2++;
//theta now equals pi/2
theta += dtheta;

	for (int i = 1; i < 2 * Iter - 1; i++){

		long double thetap = theta - dtheta;

		long double m1 = (sinl(theta) -sinl(thetap)) / (cosl(theta) - cosl(thetap));
		long double m2 = (R * sinl(theta) - rp * (long double) (2 * Iter - 1 - i) ) / (R * cosl(theta) + rp * (long double)(i) );

        long double b1 = y0 - m1 * x0;
		long double b2 = R * sinl(theta) - m2 * R * cosl(theta);

        points->x[*idx] = (b2 - b1) / (m1 - m2);
		points->y[*idx] = m1 * points->x[*idx] + y0 - m1* x0;
		xref[idx2] = points->x[*idx];
		yref[idx2] = points->y[*idx];
        x0 = points->x[*idx];
        y0 = points->y[*idx];
		(*idx)++;
		idx2++;
		theta += dtheta;
}


x0 = -rp * (long double) (2 * Iter - 1) - rp * sqrtl(2.0L);
y0 = 0.0L;
	points->x[*idx] = x0;
	points->y[*idx] = y0;
	xref[idx2] = x0;
	yref[idx2] = y0;
    
	(*idx)++;
	idx2++;
//theta now equals pi
theta += dtheta;

	for (int i = 1; i < 2 * Iter - 1; i++){

		long double thetap = theta - dtheta;

		long double m1 = (sinl(theta) -sinl(thetap)) / (cosl(theta) - cosl(thetap));
		long double m2 = (R * sinl(theta) + rp * (long double) (i) ) / (R * cosl(theta) + rp * (long double)(2 * Iter - 1 - i) );

        long double b1 = y0 - m1 * x0;
		long double b2 = R * sinl(theta) - m2 * R * cosl(theta);

        points->x[*idx] = (b2 - b1) / (m1 - m2);
		points->y[*idx] = m1 * points->x[*idx] + y0 - m1* x0;
		xref[idx2] = points->x[*idx];
		yref[idx2] = points->y[*idx];
        x0 = points->x[*idx];
        y0 = points->y[*idx];
		(*idx)++;
		idx2++;
		theta += dtheta;
}


x0 = 0.0L;
y0 = -rp * (long double) (2 * Iter - 1) - rp * sqrtl(2.0L);
	points->x[*idx] = x0;
	points->y[*idx] = y0;
	xref[idx2] = x0;
	yref[idx2] = y0;
    
	(*idx)++;
	idx2++;
//theta now equals 1.5pi
theta += dtheta;

	for (int i = 1; i < 2 * Iter - 1; i++){

		long double thetap = theta - dtheta;

		long double m1 = (sinl(theta) -sinl(thetap)) / (cosl(theta) - cosl(thetap));
		long double m2 = (R * sinl(theta) + rp * (long double) (2 * Iter - 1 - i) ) / (R * cosl(theta) - rp * (long double) (i) );

        long double b1 = y0 - m1 * x0;
		long double b2 = R * sinl(theta) - m2 * R * cosl(theta);

        points->x[*idx] = (b2 - b1) / (m1 - m2);
		points->y[*idx] = m1 * points->x[*idx] + y0 - m1* x0;
		xref[idx2] = points->x[*idx];
		yref[idx2] = points->y[*idx];
        x0 = points->x[*idx];
        y0 = points->y[*idx];
		(*idx)++;
		idx2++;
		theta += dtheta;
}
//theta is now 2pi



};

class Rocket{
	public:
};
class Hybrid: public Rocket {



};


int main(){
Rocket myrocket;
Boundary_Surface mySurface({
		.r_inner = [](double z, double phi){ return 5;},
		.r_outer = [](double z, double phi){return 6;},
});

Points mypoints;
mypoints.x.reserve(100000);
mypoints.y.reserve(100000);
mypoints.z.reserve(100000);
int idx = 0;
ButterFlyMeshPoint(5, 0.08, &idx, &mypoints);
mypoints.N = idx;
mypoints.Open();
mypoints.Write();


}
//Algorithim:
//

void nozzle_shape(long double rBig, long double rSmall, long double sharpness){

	int a[17] = {0};
	13[a] ++;
};

void pointsfill(FILE *points[3], long double r_fineness, long double rRocket, long double L, long double x0[3], void (*nozzle_shape)(long double, long double)){

	int n = 6;




}
