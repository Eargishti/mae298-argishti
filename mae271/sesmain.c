#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define l 1.0
#define m 1.0
#define g 9.8


int main(int argc, char *argv[]){
printf("");
FILE *out = fopen("lab2yi.txt", "w");
double y = 0, dt = 0.0001, tf = 10.0, t = 0.0, dtheta = 0.0;

double *theta = (double *) malloc( ((int)(tf/dt) + 1) * sizeof(double));
double *px = (double *) malloc( ((int)(tf/dt) + 1) * sizeof(double));
px[0] = 0.0;

char *endptr = NULL;

	theta[0] = 20 * M_PI / 180;
if (argc > 1){
theta[0] = strtod(argv[1], &endptr) * M_PI / 180;
}

int nfine = 100;
if (argc > 2){
	nfine = atoi(argv[2]);
}

double C = 0.1;
double w = 40.0 * M_PI;

if (argc > 3){
	C = strtod(argv[3], &endptr);
}

if (argc > 4){
	w = strtod(argv[4], &endptr) * M_PI;
}

printf("w = %.5lf\n", w/ M_PI);

double dpx = m * g * sin(theta[0]) * cos(theta[0]);


double *pth = (double *) malloc( ((int)(tf/dt) + 1) * sizeof(double));
pth[0] = 0;

/*for (int i = 0; i < (int)(tf/dt); i++) {
    double th = theta[i];

    // p_theta update first (kick)
    pth[i+1] = pth[i] + (m * g * l * sin(th)) * dt;

    // theta update (drift) using new momentum
    theta[i+1] = theta[i] + (pth[i+1] / (m * l * l)) * dt;
}*/

for (int i = 0; i < (int)(tf / dt); i++) {

    double th = theta[i];
    double ydd = -C * w * w * sin(w * t);

    double geff = g + ydd;   // = g - C w^2 sin(w t)
    double dpth = m * l * geff * sin(th);

    pth[i+1] = pth[i] + dpth * dt;

    double dtheta = pth[i+1] / (m * l * l);
    theta[i+1] = theta[i] + dtheta * dt;
    t += dt;
}




t = 0.0;

for (int i = 0; i < tf/dt; i += nfine){

double y = l * cos(theta[i]);
if (i == 0){
printf("y = %lf * %lf\n", l, cos(theta[i]));
};

fprintf(out, "%.5lf\t%.5lf\n", t, y);
t +=((double) nfine) * dt;
};

};
