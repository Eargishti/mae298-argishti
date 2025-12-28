#include <math.h>
#include <stdio.h>

int main() {

  FILE *pf = fopen("figfig", "w");

  double x1;
  double x2;
  double u1;
  double u2;
  x1 = sin(0);
  printf("1");
  double dt = 0.001;

  double x0[2] = {1.0, 0};

  double t = 0.0;
  double A_bk[2][2] = {0.0};

  A_bk[0][0] = -15.3043;
  A_bk[1][1] = -15.3043;
  A_bk[0][1] = 13.0682;
  A_bk[1][0] = 13.0682;
  x1 = 1;
  x2 = 0;

  double P1 = 0;
  double P2 = 0;
  double P3 = 0;

  double P01 = 0.6180 / 2.0;
  double P02 = 0;
  double P03 = 0.6180 / 2.0;

  P1 = P01;
  P2 = P02;
  P3 = P03;

  while (t < 10) {
    double x1_old = x1;
    x1 += A_bk[0][0] * x1 * dt + A_bk[0][1] * x2 * dt;
    x2 += A_bk[1][0] * x1_old * dt + A_bk[1][1] * x2 * dt;
    double P1_old = P1, P2_old = P2, P3_old = P3;

    P1 += 2 * P1_old * dt + 4 * P1_old * P1_old * dt +
          4 * P2_old * P2_old * dt - 1 * dt;
    P2 += 2 * P2_old * dt + 4 * P1_old * P2_old * dt + 4 * P2_old * P3_old * dt;
    P3 += 2 * P3_old * dt + 4 * P2_old * P2_old * dt +
          4 * P3_old * P3_old * dt - 1 * dt;
    if (isinf(P1))
      break;

    fprintf(pf, "%.5lf\t%.5lf\n", t, P1);
    fprintf(pf, "%.5lf\t%.5lf\n", t, P2);
    fprintf(pf, "%.5lf\t%.5lf\n", t, P3);

    t += dt;
  };
};
