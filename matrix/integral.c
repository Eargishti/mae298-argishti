#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <wchar.h>

#define BUFFER_SIZE 128
typedef struct {
  long double analytical;
  long double computed;
  long double error;

} integral;

typedef struct {
  long double weight;
  long double node;
  int i;
} Gauss;

float a4; // 4 byte

// a4 ~ 10^-7
double a5; // 8 byte

// at ~10^-15

long double a6; // 16 bytes
// a5 ~ 10^-30

long double f(long double x) {
  long double d = 1;
  long double b = 1;
  long double E = 1;
  long double L = 100;

  return ((x / L) - 1) * (x / L) / ((1 / 12.0L) * pow(2 - (x / L), 3));
}

void SetWeights(Gauss G[100], FILE *gaustxt) {
  char *endptr;
  char Buffer[BUFFER_SIZE];
  for (int i = 0; i < 50; i++) {
    fgets(Buffer, BUFFER_SIZE, gaustxt);
    G[2 * i].node = strtold(Buffer, &endptr);
    G[2 * i + 1].node = -G[2 * i].node;
    G[2 * i].weight = strtold(endptr, &endptr);
    G[2 * i + 1].weight = G[2 * i].weight;
    G[2 * i].i = 2 * i;
    G[2 * i + 1].i = 2 * i + 1;
  };
};

integral trapezoidal(long double x0, long double x1, int n) {
  integral I1;
  if (n < 0) {
    n *= -1;
  };
  I1.computed = 0;
  for (int i = 1; i <= n; i++) {

    I1.computed += f(x0 + ((long double)i - 1.0L) * (x1 - x0) / n) +
                   f(x0 + ((long double)(i)) * (x1 - x0) / n);
  };
  I1.computed *= (x1 - x0) / (2 * n);

  I1.analytical = 0;
  I1.error = 0;

  return I1;
};

long double GaussianQuadrature(long double x0, long double x1, Gauss G[100]) {
  // u = 2x/x1-x0 - 2x0 / x1 - x0 - 100

  long double result = 0.0L;
  for (int i = 0; i < 100; i++) {
    result += ((x1 - x0) / 2.0L) * G[i].weight *
              f((G[i].node + 1) * (x1 - x0) / 2 + x0);
  };

  return result;
};

int main() {
  long double x0, x1;
  int n;
  integral I1;
  x0 = 0;
  x1 = 100;
  n = 15;
  printf("f(x) = ((x/L) - 1) * (x/L) / ( (1/12.0L) * pow(2 - (x/L), 3))\n");
  I1 = trapezoidal(x0, x1, n);
  printf("x0 = %.5Lf, x1 = %.5Lf, n = %d\n", x0, x1, n);
  printf("I = %.25Lf\n", I1.computed);
  FILE *f1;
  f1 = fopen("ggauspoints", "r");
  Gauss G1[100];
  SetWeights(G1, f1);
  long double I0;

  I0 = GaussianQuadrature(x0, x1, G1);
  printf("x0 = %.5Lf, x1 = %.5Lf, n = %d\n", x0, x1, n);
  printf("I = %.25Lf\n", I0);
};

/*wc = fgetwc(f1);

fprintf(stdout, "plus or minus from file = %04X, number = %d\n", (unsigned)wc,
        (unsigned)wc);
fprintf(stdout, "character = %lc\n", wc); */
