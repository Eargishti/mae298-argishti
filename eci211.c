#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
  long double area[3] = {15000, 20000, 10000};

  long double length[3] = {8000, 5657, 10928};
  long double angle[3] = {210, 315, 0.0};
  int neIem = 3;
  long double klocal[4][4] = {0};
  const long double E = 200;

  klocal[0][0] = 1;
  klocal[2][0] = -1;
  klocal[0][2] = -1;
  klocal[2][2] = 1;
  long double tran[4][4] = {0};

  long double EAL;
  long double theta;
  long double C[3];
  long double S[3];

  for (int i = 0; i < 4; i++) {

    EAL = E * area[i] / length[i];
    theta = angle[i] * M_PI / 180;
    C[i] = cosl(theta);
    S[i] = sinl(theta);

    for (int j = 0; j < 4; j++) {
      klocal[i][j] *= EAL;
    };
  };

  for (int i = 0; i < neIem; i++) {
  };
};
