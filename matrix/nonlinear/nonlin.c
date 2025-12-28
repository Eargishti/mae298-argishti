#include <math.h>
#include <stdio.h>

double fx(double uy, double ux) {
  double E = 29000;
  double Acb = 0.4;
  double Aab = 3.1415926353589 * (1.0 / 16.0) * (1 / 16.0);
  double fff = (-E * Acb * (sqrt(pow(60 + uy, 2) + ux * ux) - 60) * ux /
                (pow(60 + uy, 2) + ux * ux)) -
               E * Aab *
                   (sqrt(pow(60 + uy, 2) + pow(60 + ux, 2)) - 60 * sqrt(2)) *
                   (60 + ux) / (pow(60 + ux, 2) + pow(60 + uy, 2));
  printf("\n");
  return fff;
};

double fy(double uy, double ux) {
  double E = 29000;
  double Acb = 0.4;
  double Aab = 3.1415926353589 * (1.0 / 16.0) * (1 / 16.0);
  double fff = (-E * Acb * (sqrt(pow(60 + uy, 2) + ux * ux) - 60) * (60 + uy) /
                (pow(60 + uy, 2) + ux * ux)) -
               E * Aab *
                   (sqrt(pow(60 + uy, 2) + pow(60 + ux, 2)) - 60 * sqrt(2)) *
                   (60 + uy) / (pow(60 + ux, 2) + pow(60 + uy, 2));
  printf("\n");
  return fff;
};

int main() {
  double E = 29000;
  double Acb = 0.4;
  double Aab = 3.1415926353589 * (1.0 / 16.0) * (1 / 16.0);
  double ux = 0.06911717373;

  ux = 0.088063876712;
  double uy = -0.077436982005;
  uy = -0.077664385497;
  printf("(ux,uy) = (%.12lf, %.12lf)\n", ux, uy);
  printf("(fx,fy) = (%.9lf, %.9lf)\n", fx(uy, ux), fy(uy, ux));
  printf("Axial force_ab = %.9lf\n",
         E * Aab * (sqrt(pow(60 + ux, 2) + pow(60 + uy, 2)) - 60 * sqrt(2)) /
             sqrt(pow(60 + ux, 2) + pow(60 + uy, 2)));
  // printf("Axial force_bc = %.9lf\n", E * Acb *(sqrt(pow(60+ux,2) +
  // pow(60+uy,2)) - 60 * sqrt(2)) /  sqrt(pow(60+ux,2) + pow(60+uy,2)) );
};
// y=-0.077664385497
// x=0.088063876712
