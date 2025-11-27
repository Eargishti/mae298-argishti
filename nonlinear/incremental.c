#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const double del_uy = 0.001;
const double E = 29000;
const double L1 = 16 * 12;
const double L2 = 48;
const double A = 6;

double fbar(double A, double L, double uy) {
  double ubar = sqrt(L1 * L1 + (L2 + uy) * (L2 + uy)) - sqrt(L1 * L1 + L2 * L2);
  return E * ubar * A / L;
};
double fy(double A, double L, double uy) {

  return fbar(A, L, uy) * (L2 + uy) / sqrt(L1 * L1 + (L2 + uy) * (L2 + uy));
};

int main() {
  double L0 = sqrt(L1 * L1 + L2 * L2);
  double uy = 0;
  double Fy = 0;
  double error = 0.0001;
  double diff = 150;
  double Fbar = 0;

  while (diff > error) {
    uy = uy + del_uy;
    double Lstar = sqrt(L1 * L1 + (L2 + uy) * (L2 + uy));
    Fbar = fbar(A, Lstar, uy);
    Fy = fy(A, Lstar, uy);
    diff = 150 - Fy;
  };
  printf("\nuy = %f\n", uy);
};
