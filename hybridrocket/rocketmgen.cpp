#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <functional>
#include "headerg.cpp"


std::string Edgar = "hi";

struct Radial_Functions {

std::function<double(double zeta, double phi)> r_inner;
std::function<double(double zeta, double phi)> r_outer;
double z_length = 0;

};


struct Circular_let{
std::function<double(double phi)> r_inner;
double z;
};




class Boundary_Surface{
private:
Radial_Functions R1;
Circular_let C1;
std::function<double(double zeta)>  phi_min;
std::function<double(double zeta)>  phi_max;
long double z_min, z_max;

public:


enum shape{
CIRCLULAR_PLATE,
CYLINDRICAL_SURFACE,
} Shape;

Boundary_Surface(Radial_Functions p) : R1(p){Shape = CYLINDRICAL_SURFACE;}

void set_function(std::function<double(double, double)> Boundary_Surface::*target, std::function<double(double, double)> source){
this->*target = source;
};

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
		.phi_min = [](double z){ return 0; },
		.phi_max = [](double z){ return M_PI; }
});

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
