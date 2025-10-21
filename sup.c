#include <math.h>
#include <stdio.h>

typedef struct {
  long double analytical;
  long double computed;
  long double error;

} integral;

long double f(long double x) {
  long double d = 1;
  long double b = 1;
  long double E = 1;
  long double L = 100;

  return ((x / L) - 1) * (x / L) / ((1 / 12.0L) * pow(2 - (x / L), 3));
}

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

int main() {
  long double x0, x1;
  int n;
  integral I1;
  x0 = 0;
  x1 = 100;
  n = 15000;
  printf("f(x) = ((x/L) - 1) * (x/L) / ( (1/12.0L) * pow(2 - (x/L), 3))\n");
  I1 = trapezoidal(x0, x1, n);
  printf("x0 = %.5Lf, x1 = %.5Lf, n = %d\n", x0, x1, n);
  printf("I = %.9Lf\n", I1.computed);
};
