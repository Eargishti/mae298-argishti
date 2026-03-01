#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#define N 4E+06
#define pi M_PI
#define points 2000

double randn() {
    // Box-Muller transform
    double u1 = (rand() + 1.0) / (RAND_MAX + 2.0);
    double u2 = (rand() + 1.0) / (RAND_MAX + 2.0);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}


double interp1(double x, const double *y, int n_pts, double delta_x)
{
    if (x <= 0.0) return y[0];
    if (x >= n_pts * delta_x) return y[n_pts];

    int i = (int)(x / delta_x);

    double x_i = i * delta_x;
    double t = (x - x_i) / delta_x;

    return (1.0 - t) * y[i] + t * y[i + 1];
}


int main(int argc, char *argv[]){
    srand((unsigned)time(NULL));
	char *endptr;
double zeta_s = 0.7;
double zeta_c = 0.7;
if (argc != 1){
zeta_s = strtod(argv[1], &endptr);
zeta_c = strtod(argv[2], &endptr);
};
double delta_x = 0.5;
double Length = 500;
double U = 40  *  .46;
double g = 9.8;
double mtot = 3000 / 2.2;
double msmus = 5;
double mus = mtot / (1 + msmus);
double ms = mtot - mus;
double ws = 2 * pi * 1.2;
double wwh = 2 * pi * 8;
double Rw = .001;
double T = 10;
double ma = .03 * ms;
double wa = 2 * pi * 15;

double kt = mus * wwh  * wwh;
double ka = ma *  wa  *  wa;
double ks = ms * ws  *  ws;
double ba = 2 * .1 * wa * ma;
double bs = 2 * zeta_s * ws * ms;
double bc = 2 * zeta_c * ws * ms;




int n_pts = (int) (Length / delta_x);
    double *X_i      = malloc((n_pts + 1) * sizeof(double));
    double *slope_raw = malloc((n_pts + 1) * sizeof(double));
    double *slope_i   = malloc((n_pts + 1) * sizeof(double));

  for (int i = 0; i <= n_pts; i++) {
        X_i[i] = i * delta_x;
    }

    double mean = 0.0;
    for (int i = 0; i <= n_pts; i++) {
        slope_raw[i] = randn();
        mean += slope_raw[i];
    }
    mean /= (n_pts + 1);

    for (int i = 0; i <= n_pts; i++) {
        slope_i[i] = 0.005 * (slope_raw[i] - mean);
    }

double dpus = 0, dps = 0, dpa = 0, dqt = 0, dqs = 0, dqa = 0;
double X = 0;
double slope = 0;
//X = U *t
//slope = interpl(X_i, slope_i, X);
//v_i = U * slope;

double vi = 0;

double t = 0;
double tfinal = 29;
double *ps = (double *)  calloc(N , sizeof(double));
double *pus = (double *) calloc(N , sizeof(double));
double *pa = (double *)  calloc(N , sizeof(double));


double *qt = (double *) calloc(N , sizeof(double));
double *qs = (double *) calloc(N , sizeof(double));
double *qa = (double *) calloc(N , sizeof(double));
double *ic = (double *) calloc(N , sizeof(double));

double dt = tfinal / N;
int i = 0;

while (t < tfinal){

	ic[i] = ps[i] * bc /(T * ms);
	X = U * t;
	slope = interp1(X, slope_i, n_pts, delta_x);
	vi = U * slope;

dpus = -mus * g + qt[i] * kt - qs[i] * ks - (pus[i] / mus - ps[i] /ms) * bs;
dps = - ms * g + qs[i] * ks + (pus[i] / mus - ps[i] / ms ) * bs - qa[i] * ka - (ps[i] / ms - pa[i] / ma) * ba - T * ic[i];
dpa = -ma * g + qa[i] * ka + (ps[i] / ms - pa[i] /ma) * ba + T * ic[i];
dqt = vi - pus[i] / mus;
dqs = pus[i] / mus - ps[i] / ms;
dqa = ps[i]/ms - pa[i]/ma;

 ps[i + 1] = dps * dt + ps[i];
pus[i + 1] = dpus * dt + pus[i];
 pa[i + 1] = dpa * dt + pa[i];
 qt[i + 1] = dqt * dt + qt[i];
 qs[i + 1] = dqs * dt + qs[i];
 qa[i + 1] = dqa * dt + qa[i];


t += dt;
i ++;
};

t = 0;
int ror = N / points;
FILE *out = fopen("as_vs_t","w"); 
FILE *out2 = fopen("qa_vs_t","w"); 
FILE *out3 = fopen("P_vs_t","w"); 
FILE *out4 = fopen("aus_vs_t","w"); 




for (i = 0; i < points ; i++){
fprintf(out, "%.6lf\t%.6lf\n", t, ps[i*ror] / ms);
fprintf(out2, "%.6lf\t%.6lf\n", t, qa[i*ror] / ms);
fprintf(out3, "%.6lf\t%.6lf\n", t, ic[i * ror] * ic[i * ror] * Rw);
fprintf(out4, "%.6lf\t%.6lf\n", t, pus[i*ror] / mus);
t += tfinal / points;
};
double P_rms = 0;
double as_rms = 0;

for (i = 0; i < N; i++){
P_rms += pow(ic[i] * ic[i] * Rw, 2.0) ;
as_rms += pow(ps[i] / ms, 2.0);
};

P_rms = sqrt(P_rms);
as_rms = sqrt(as_rms);
P_rms /= sqrt(N);
as_rms /= sqrt(N);
fprintf(stdout, "P_rms = %.6lf\nas_rms = %.6lf\n", P_rms, as_rms);



};
