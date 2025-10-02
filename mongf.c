#include <locale.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFFER_SIZE 128
/*long double a[] = {0.2969, -.1260, -.3516, .2843, -.1015};
long double m = 0.02;
long double p = 0.4;
long double T = 0.11;
long double kk = m / (p * p);
long double jj = m / ((1 - p) * (1 - p));
long double cpar = 0.0041;
long double t4 = 1.0089304113; */

const long double delz = 0.05L;
const long double tol = 1e-15L;
typedef struct {
  long double m;
  long double T;
  long double p;
  long double kk;
  long double jj;
  long double a[5];
  long double c;
} aerofoil;

typedef struct {
  long double xL;
  long double xR;
  long double yU;
  long double yB;
  int mp;
  int r_num;
  int r_dem;
  int gp;
  long double Cpar;
  int Tpar;

} MeshFineness;

long double f1(long double t, aerofoil a) {
  return (a.kk * (2 * a.p * (t / a.c) - ((t / a.c) * (t / a.c))));
}

long double f2(long double t, aerofoil a) {
  return (a.jj * (1 - 2 * a.p + 2 * a.p * (t / a.c) - (t / a.c) * (t / a.c)));
}

long double f3(long double t, aerofoil a) {
  return (5 * a.T *
          (a.a[0] * sqrt(t / a.c) + a.a[1] * (t / a.c) +
           a.a[2] * pow((t / a.c), 2.0) + a.a[3] * pow((t / a.c), 3.0) +
           a.a[4] * pow((t / a.c), 4.0)));
}

long double theta1(long double t, aerofoil a) {
  return (atan(2 * a.kk * (a.p / a.c - t / (a.c * a.c))));
}

long double theta2(long double t, aerofoil a) {
  return (atan(2 * a.jj * (a.p / a.c - t / (a.c * a.c))));
}

long double f1prime(long double t, aerofoil a) {
  return (2 * a.kk * (a.p / a.c - t / (a.c * a.c)));
}

long double f2prime(long double t, aerofoil a) {
  return (2 * a.jj * (a.p / a.c - t / (a.c * a.c)));
}

long double f3prime(long double t, aerofoil a) {
  return (5 * a.T *
          (a.a[0] * (1 / (2 * sqrt(a.c * t))) + a.a[1] / a.c +
           2 * a.a[2] * t / (a.c * a.c) + 3 * a.a[3] * t * t / (pow(a.c, 3.0)) +
           4 * a.a[4] * pow(t, 3.0) / (pow(a.c, 4.0))));
}

long double theta1prime(long double t, aerofoil a) {
  return (-2 * a.kk /
          ((a.c * a.c) *
           (1 + pow((2 * a.kk * (a.p / a.c - t / (a.c * a.c))), 2.0))));
}

long double theta2prime(long double t, aerofoil a) {
  return (-2 * a.jj /
          ((a.c * a.c) *
           (1 + pow((2 * a.jj * (a.p / a.c - t / (a.c * a.c))), 2.0))));
}

long double dydx1(long double t, aerofoil a) {
  return ((f1prime(t, a) + f3prime(t, a) * cos(theta1(t, a)) -
           theta1prime(t, a) * f3(t, a) * sin(theta1(t, a))) /
          (1.0 - f3prime(t, a) * sin(theta1(t, a)) -
           theta1prime(t, a) * f3(t, a) * cos(theta1(t, a))));
}

long double dydx2(long double t, aerofoil a) {
  return (f1prime(t, a) - f3prime(t, a) * cos(theta1(t, a)) +
          theta1prime(t, a) * f3(t, a) * sin(theta1(t, a))) /
         (1 + f3prime(t, a) * sin(theta1(t, a)) +
          theta1prime(t, a) * f3(t, a) * cos(theta1(t, a)));
}

long double dydx3(long double t, aerofoil a) {
  return (f2prime(t, a) + f3prime(t, a) * cos(theta2(t, a)) -
          theta2prime(t, a) * f3(t, a) * sin(theta2(t, a))) /
         (1 - f3prime(t, a) * sin(theta2(t, a)) -
          theta2prime(t, a) * f3(t, a) * cos(theta2(t, a)));
}

long double dydx4(long double t, aerofoil a) {
  return (f2prime(t, a) - f3prime(t, a) * cos(theta2(t, a)) +
          theta2prime(t, a) * f3(t, a) * sin(theta2(t, a))) /
         (1 + f3prime(t, a) * sin(theta2(t, a)) +
          theta2prime(t, a) * f3(t, a) * cos(theta2(t, a)));
}

// Root functions. Zeroroot is t0, Firstroot is t1, Secondroot is t2
long double Zeroroot(long double testval, aerofoil aero) {
  return (1 - f3prime(testval, aero) * sin(theta1(testval, aero)) -
          f3(testval, aero) * theta1prime(testval, aero) *
              cos(theta1(testval, aero))) /
         (f1prime(testval, aero) +
          f3prime(testval, aero) * cos(theta1(testval, aero)) -
          theta1prime(testval, aero) * f3(testval, aero) *
              sin(theta1(testval, aero)));
};

long double Firstroot(long double t, aerofoil a, long double xL,
                      long double yU) {
  // f1(t) + f3(t) * cos (theta1(t) );
  // t - f3(t) * sin (theta1(t) );

  return ((yU - (f1(t, a) + f3(t, a) * cos(theta1(t, a)))) /
          (xL - t + f3(t, a) * sin(theta1(t, a)))) +
         Zeroroot(t, a);
};
long double Secondroot(long double t, aerofoil a, long double xL,
                       long double yB) {

  return ((yB - (f1(t, a) - f3(t, a) * cos(theta1(t, a)))) /
          (xL - t - f3(t, a) * sin(theta1(t, a)))) +
         (1 + f3prime(t, a) * sin(theta1(t, a)) +
          f3(t, a) * theta1prime(t, a) * cos(theta1(t, a))) /
             (f1prime(t, a) - f3prime(t, a) * cos(theta1(t, a)) +
              theta1prime(t, a) * f3(t, a) * sin(theta1(t, a)));
};

long double Thirdroot(long double t, aerofoil a) {

  return 2 * f3(t, a) * cos(theta2(t, a));
};

long double yroot(long double t, aerofoil a, long double yval) {

  return (f1(t, a) + f3(t, a) * cosl(theta1(t, a)) - yval);
};

typedef struct foilpoint {
  long double *x;
  long double *y;
  long double *slope;
  long double *associatedT;
  int *label;
  int *section;
  int count;

} foilpoint;

typedef struct trace {
  long double *x;
  long double y;
  int *label;
  int size;

} trace;

typedef struct Tarray {
  int size;

  long double *T;

} Tarray;

typedef struct cell {
  int label;
  int faces[6];

} cell;

void GfaceFILL(FILE *faces[3], int o4, int o3, int o2, int o1, int *flabel,
               int bcode, cell *oker, int faceN, int *clabel) {
  faceN++;
  oker->faces[faceN] = *flabel;
  oker->label = *clabel;
  o1++;
  o2++;
  o3++;
  o4++;
  fprintf(faces[0], "4(%d %d %d %d)\n", o1, o2, o3, o4);
  // fprintf(faces[0], "( %d %8d %8d %8d )\n", o1, o2, o3, o4);
  fprintf(faces[1], "( %d %8d %8d %8d )\t%8d\t%8d\n", o1, o2, o3, o4, bcode,
          *flabel);
  // fprintf(faces[2], "%d\t%d\t%d\t%d\t%8d\t%8d\n", o1, o2, o3, o4, bcode,
  // *flabel);
  fwrite(&o1, sizeof(int), 1, faces[2]);
  fwrite(&o2, sizeof(int), 1, faces[2]);
  fwrite(&o3, sizeof(int), 1, faces[2]);
  fwrite(&o4, sizeof(int), 1, faces[2]);
  fwrite(&bcode, sizeof(int), 1, faces[2]);
  fwrite(flabel, sizeof(int), 1, faces[2]);

  (*flabel)++;

  (*clabel)++;
};
typedef enum {
  AEROFOIL,
  NORTH,
  SOUTH,
  EAST,
  WEST,
  IN,
  OUT,
  INLET,
  OUTLET,
  TOP,
  BOTTOM,
  REVERSE,
  TOP2,
  BOTTOM2
} facecode;

typedef enum {
  FrontAndBack,
  Aerofoil,
  Inlet,
  Outlet,
  Top1,
  Bottom1,
  Top2,
  Bottom2
} BoundaryPatch;
//    			  0		   1	  2    		  3         4 5
//    6      7

void ownFILL(FILE *owner[2], FILE *neighbour[2], int clabelo, int clabeln,
             facecode *bc) {

  if (*bc != AEROFOIL && *bc != IN && *bc != OUT && *bc != INLET &&
      *bc != OUTLET && *bc != TOP && *bc != BOTTOM && *bc != TOP2 &&
      *bc != BOTTOM2) {
    fprintf(neighbour[0], "%d\n", clabeln);
    fwrite(&clabeln, sizeof(int), 1, neighbour[1]);
    fwrite(bc, sizeof(int), 1, neighbour[1]);
  };

  fprintf(owner[0], "%d\n", clabelo);
  // fprintf(owner[1], "%d\t%d\n", clabelo, *bc);
  fwrite(&clabelo, sizeof(int), 1, owner[1]);
  fwrite(bc, sizeof(int), 1, owner[1]);
};

void faceFILL(FILE *faces[3], int o1, int o2, int o3, int o4, int *flabel,
              int bcode, cell *oker, int faceN, int *clabel) {
  oker->faces[faceN] = *flabel;
  oker->label = *clabel;

  fprintf(faces[0], "4(%d %d %d %d)\n", o1, o2, o3, o4);
  // fprintf(faces[0], "( %d %8d %8d %8d )\n", o1, o2, o3, o4);
  fprintf(faces[1], "( %d %8d %8d %8d )\t%8d\t%8d\n", o1, o2, o3, o4, bcode,
          *flabel);
  // fprintf(faces[2], "%d\t%d\t%d\t%d\t%8d\t%8d\n", o1, o2, o3, o4, bcode,
  // *flabel);
  fwrite(&o1, sizeof(int), 1, faces[2]);
  fwrite(&o2, sizeof(int), 1, faces[2]);
  fwrite(&o3, sizeof(int), 1, faces[2]);
  fwrite(&o4, sizeof(int), 1, faces[2]);
  fwrite(&bcode, sizeof(int), 1, faces[2]);
  fwrite(flabel, sizeof(int), 1, faces[2]);

  (*flabel)++;

  (*clabel)++;
  //  cells[0] = fopen("cells", "w");
  //  cells[1] = fopen("diagnosis", "w");
  //  cells[2] = fopen("cellcheck", "w");
};

void onlyFILL(FILE *faces[3], int o1, int o2, int o3, int o4, int *flabel,
              int bcode) {
  fprintf(faces[0], "4(%d %d %d %d)\n", o1, o2, o3, o4);
  // fprintf(faces[0], "( %d %8d %8d %8d )\n", o1, o2, o3, o4);
  fprintf(faces[1], "( %d %8d %8d %8d )\t%8d\t%8d\n", o1, o2, o3, o4, bcode,
          *flabel);
  // fprintf(faces[2], "%d\t%d\t%d\t%d\t%8d\t%8d\n", o1, o2, o3, o4, bcode,
  // *flabel);
  fwrite(&o1, sizeof(int), 1, faces[2]);
  fwrite(&o2, sizeof(int), 1, faces[2]);
  fwrite(&o3, sizeof(int), 1, faces[2]);
  fwrite(&o4, sizeof(int), 1, faces[2]);
  fwrite(&bcode, sizeof(int), 1, faces[2]);
  fwrite(flabel, sizeof(int), 1, faces[2]);
};

void MaxFILL(FILE *f1, int count) {
  char buff1[128];
  for (int i = 0; i < 15; i++) {
    fgets(buff1, sizeof(buff1), f1);
  };

  fprintf(f1, "%d\n", count);
};
void headsF(FILE *f1, const char *sting, const char *stong) {

  fprintf(f1, "");
  fprintf(f1, " /*--------------------------------*- C++ "
              "-*----------------------------------*\\\n");
  fprintf(f1, "  =========                 |\n");

  fprintf(f1, "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD "
              "Toolbox\n");

  fprintf(f1,
          "   \\\\    /   O peration     | Website:  https://openfoam.org\n");
  fprintf(f1, "    \\\\  /    A nd           | Version:  dev\n");
  fprintf(f1, "     \\\\/     M anipulation  |\n");
  fprintf(f1, "\\*-------------------------------------------------------------"
              "--------------*/\n");
  fprintf(f1, "FoamFile\n");
  fprintf(f1, "{\n");
  fprintf(f1, "    format      ascii;\n");
  fprintf(f1, "    class       %s;\n", sting);
  fprintf(f1, "    object      %s;\n", stong);
  fprintf(f1, "}\n");
  fprintf(f1, " // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * "
              "* * * * * * * //");
};

void headerfill() {
  FILE *points = fopen("points", "w");

  FILE *faces = fopen("faces", "w");
  FILE *cells = fopen("cells", "w");
  FILE *owner = fopen("owner", "w");
  FILE *neighbour = fopen("neighbour", "w");
  FILE *boundary = fopen("boundary", "w");
  headsF(points, "vectorField", "points");
  headsF(faces, "faceList", "faces");
  headsF(cells, "cellList", "cells");
  headsF(owner, "labelList", "owner");
  headsF(neighbour, "labelList", "neighbour");
  headsF(boundary, "polyBoundaryMesh", "boundary");
  fclose(points);
  fclose(faces);
  fclose(cells);
  fclose(owner);
  fclose(neighbour);
  fclose(boundary);
};

typedef enum { four, TopAndBottom } OPTIONB;
void WriteToPoint(const void *restrict x, const void *restrict y,
                  const void *restrict z0, const void *restrict label,
                  FILE *pointoutput) {

  fwrite(x, sizeof(long double), 1, pointoutput);
  fwrite(y, sizeof(long double), 1, pointoutput);
  fwrite(z0, sizeof(long double), 1, pointoutput);
  fwrite(label, sizeof(int), 1, pointoutput);
};

// typedef enum { FrontAndBack, Aerofoil, Inlet, Outlet, Top1, Bottom1, Top2,
// Bottom2 } BoundaryPatch;
const char *GetBoundaryName(BoundaryPatch Patch1) {
  switch (Patch1) {
  case FrontAndBack:
    return "FrontAndBack";
  case Aerofoil:
    return "Aerofoil";
  case Inlet:
    return "Inlet";
  case Outlet:
    return "Outlet";
  case Top1:
    return "Top1";
  case Bottom1:
    return "Bottom1";
  case Top2:
    return "Top2";
  case Bottom2:
    return "Bottom2";
  }
};

void linestack(FILE *pointoutput[4], foilpoint *foil, foilpoint *zoinks,
               aerofoil a, long double cpar, Tarray array, int ncrit,
               long double xL, int mp, int r_num, int r_dem, long double xR,
               long double yU, long double yB, int backtrace, trace *trace1,
               long double slopeshift, FILE *faces[3], FILE *cells[3],
               FILE *owner[2], FILE *neighbour[2], FILE *boundar[2],
               OPTIONB opt1) {
  pointoutput[0] = fopen("figaro1.txt", "a");
  pointoutput[1] = fopen("figarosafe.txt", "a");
  pointoutput[2] = fopen("points", "a");
  pointoutput[3] = fopen("pointcheck", "a");
  if ((pointoutput[0] == NULL || pointoutput[1] == NULL)) {
    printf("\nERROR OPENING FILE\n");
  };

  long double x1, y1, x2, deltax;
  x2 = 0;
  x1 = 0.0;
  y1 = 0.0;
  int n = 0;
  int label1 = foil->count + backtrace;
  long double ytail = foil->y[array.size - 1];
  printf("n = %d\n", n);
  printf("foil->count\n = %d", foil->count);
  printf("label1\n = %d", label1);
  int ntip = array.size - 1;
  printf("Within linestack, ntip = %d\n", ntip);
  printf("Associated T at ntip = %.15Lf\n", foil->associatedT[ntip]);
  int j = 1;
  long double cpari = cpar;
  int tracecount = 1;

  if (ncrit == 1) {
    ncrit = 1;
  };
  printf("\nfoil count = %d", foil->count);
  printf("\n ncrit = %d", ncrit);
  printf("\n xtail = %.15Lf", foil->x[array.size - 1]);
  printf("\n ntip = %d", ntip);
  printf("\n foil label tail = %d", foil->label[ntip]);
  printf("\n t at ncrit = %.15Lf", foil->associatedT[ncrit]);
  printf("\n t val of corner = %.15Lf", foil->associatedT[ncrit + mp]);
  printf("\n t val of bottom left corner = %.15Lf",
         foil->associatedT[ncrit + (int)(2 * mp / 3)]);
  printf("\n last element of foil count is %.15Lf",
         foil->associatedT[foil->count - ncrit - (int)(2 * mp / 3)]);

  while (x2 > (3.0 / 5.0) * xL) {

    while (n < foil->count) {
      if (foil->section[n] == 1) {
        if (foil->label[n] <= ncrit || foil->slope[n] < 0.0) {
          foil->x[n] = foil->x[n] - cpar * cosl(atanl(foil->slope[n]));
          foil->y[n] -= cpar * sinl(atanl(foil->slope[n]));
        } else {
          foil->x[n] += cpar * cosl(atanl(foil->slope[n]));
          foil->y[n] += cpar * sinl(atanl(foil->slope[n]));
        }
        fprintf(pointoutput[0], "%.15Lf \t %.15Lf\n", foil->x[n], foil->y[n]);

        fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[n],
                foil->y[n], 0.0L, 2 * label1);
        fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[n],
                foil->y[n], delz, 2 * label1 + 1);

        WriteToPoint(&foil->x[n], &foil->y[n], &(long double){0.0L},
                     &(int){2 * label1}, pointoutput[3]);
        WriteToPoint(&foil->x[n], &foil->y[n], &(long double){delz},
                     &(int){2 * label1 + 1}, pointoutput[3]);

        fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[n],
                foil->y[n], 0.0L);
        fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[n],
                foil->y[n], delz);

        if (foil->label[n] == ncrit) {
          x2 = foil->x[n];
          // printf("We're stuck here %.15Lf > %.15Lf slope = %.15Lf dx = %.15Lf
          // \n", x2, (3.0 / 5.0) * xL, foil->slope[n], cpar *
          // cosl(atanl(foil->slope[n])));
        };
      };
      if (foil->section[n] == 3) {
        if (foil->slope[n] > 0.0 && (foil->label[n] != ntip)) {
          foil->x[n] += cpar * cosl(atanl(foil->slope[n]));
          foil->y[n] += cpar * sinl(atanl(foil->slope[n]));
          fprintf(pointoutput[0], "%.15Lf \t %.15Lf\n", foil->x[n], foil->y[n]);

          fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[n],
                  foil->y[n], 0.0L, 2 * label1);
          fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[n],
                  foil->y[n], delz, 2 * label1 + 1);

          WriteToPoint(&foil->x[n], &foil->y[n], &(long double){0.0L},
                       &(int){2 * label1}, pointoutput[3]);
          WriteToPoint(&foil->x[n], &foil->y[n], &(long double){delz},
                       &(int){2 * label1 + 1}, pointoutput[3]);

          fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[n],
                  foil->y[n], 0.0L);
          fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[n],
                  foil->y[n], delz);
        };
        if (foil->slope[n] < 0.0 && (foil->label[n] != ntip)) {
          foil->x[n] -= cpar * cosl(atanl(foil->slope[n]));
          foil->y[n] -= cpar * sinl(atanl(foil->slope[n]));
          fprintf(pointoutput[0], "%.15Lf \t %.15Lf\n", foil->x[n], foil->y[n]);

          fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[n],
                  foil->y[n], 0.0L, 2 * label1);
          fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[n],
                  foil->y[n], delz, 2 * label1 + 1);

          WriteToPoint(&foil->x[n], &foil->y[n], &(long double){0.0L},
                       &(int){2 * label1}, pointoutput[3]);
          WriteToPoint(&foil->x[n], &foil->y[n], &(long double){delz},
                       &(int){2 * label1 + 1}, pointoutput[3]);

          fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[n],
                  foil->y[n], 0.0L);
          fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[n],
                  foil->y[n], delz);
        };
        if (foil->label[n] == ntip) {
          foil->x[n] += cpar * cosl(atanl(foil->slope[n]));
          foil->y[n] += cpar * sinl(atanl(foil->slope[n]));

          fprintf(pointoutput[0], "%.15Lf \t %.15Lf\n", foil->x[n], foil->y[n]);

          fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[n],
                  foil->y[n], 0.0L, 2 * label1);
          fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[n],
                  foil->y[n], delz, 2 * label1 + 1);

          WriteToPoint(&foil->x[n], &foil->y[n], &(long double){0.0L},
                       &(int){2 * label1}, pointoutput[3]);
          WriteToPoint(&foil->x[n], &foil->y[n], &(long double){delz},
                       &(int){2 * label1 + 1}, pointoutput[3]);

          fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[n],
                  foil->y[n], 0.0L);
          fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[n],
                  foil->y[n], delz);

          label1++;
          for (int i = 1; i <= backtrace; i++) {
            long double deltax = (xR - foil->x[n]) / ((long double)(backtrace));
            trace1->x[i - 1] = ((long double)(i)) * deltax + foil->x[n];
            trace1->y = foil->y[n];
            trace1->label[i - 1] = label1;
            fprintf(pointoutput[0], "%.15Lf\t%.15Lf\n", trace1->x[i - 1],
                    trace1->y);
            fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                    trace1->x[i - 1], trace1->y, 0.0L, 2 * label1);
            fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                    trace1->x[i - 1], trace1->y, delz, 2 * label1 + 1);

            WriteToPoint(&trace1->x[i - 1], &trace1->y, &(long double){0.0L},
                         &(int){2 * label1}, pointoutput[3]);
            WriteToPoint(&trace1->x[i - 1], &trace1->y, &(long double){delz},
                         &(int){2 * label1 + 1}, pointoutput[3]);

            fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n",
                    trace1->x[i - 1], trace1->y, 0.0L);
            fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n",
                    trace1->x[i - 1], trace1->y, delz);
            label1++;
          };

          tracecount++;
          ytail -= cpar * sinl(atanl(foil->slope[n]));
          fprintf(pointoutput[0], "%.15Lf \t %.15Lf\n", foil->x[n], ytail);

          fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[n],
                  ytail, 0.0L, 2 * label1);
          fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[n],
                  ytail, delz, 2 * label1 + 1);

          WriteToPoint(&foil->x[n], &ytail, &(long double){0.0L},
                       &(int){2 * label1}, pointoutput[3]);
          WriteToPoint(&foil->x[n], &ytail, &(long double){delz},
                       &(int){2 * label1 + 1}, pointoutput[3]);

          fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[n], ytail,
                  0.0L);
          fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[n], ytail,
                  delz);

          label1++;
          for (int i = 1; i <= backtrace; i++) {
            long double deltax = (xR - foil->x[n]) / ((long double)(backtrace));
            trace1->x[i - 1] = ((long double)(i)) * deltax + foil->x[n];
            trace1->y = ytail;
            trace1->label[i - 1] = label1;
            fprintf(pointoutput[0], "%.15Lf\t%.15Lf\n", trace1->x[i - 1],
                    trace1->y);

            fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                    trace1->x[i - 1], trace1->y, 0.0L, 2 * label1);
            fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                    trace1->x[i - 1], trace1->y, delz, 2 * label1 + 1);

            WriteToPoint(&trace1->x[i - 1], &trace1->y, &(long double){0.0L},
                         &(int){2 * label1}, pointoutput[3]);
            WriteToPoint(&trace1->x[i - 1], &trace1->y, &(long double){delz},
                         &(int){2 * label1 + 1}, pointoutput[3]);

            fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n",
                    trace1->x[i - 1], trace1->y, 0.0L);
            fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n",
                    trace1->x[i - 1], trace1->y, delz);

            label1++;
          };
          label1--;

          tracecount++;
        };
      };
      if (foil->section[n] == 4) {

        long double p1 = (foil->associatedT[n] - a.c * a.p) /
                         (foil->associatedT[ntip] - a.c * a.p);
        long double p2 = 1.0 - p1;
        long double tempslope =
            (0.0 - foil->slope[ntip] - slopeshift) * p1 + p2 * foil->slope[n];
        if (tempslope > 0.0) {
          foil->x[n] -= cpar * cosl(atanl(tempslope));
          foil->y[n] -= cpar * sinl(atanl(tempslope));
          fprintf(pointoutput[0], "%.15Lf \t %.15Lf\n", foil->x[n], foil->y[n]);

          fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[n],
                  foil->y[n], 0.0L, 2 * label1);
          fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[n],
                  foil->y[n], delz, 2 * label1 + 1);

          WriteToPoint(&foil->x[n], &foil->y[n], &(long double){0.0L},
                       &(int){2 * label1}, pointoutput[3]);
          WriteToPoint(&foil->x[n], &foil->y[n], &(long double){delz},
                       &(int){2 * label1 + 1}, pointoutput[3]);

          fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[n],
                  foil->y[n], 0.0L);
          fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[n],
                  foil->y[n], delz);

        } else {
          foil->x[n] += cpar * cosl(atanl(tempslope));
          foil->y[n] += cpar * sinl(atanl(tempslope));
          fprintf(pointoutput[0], "%.15Lf \t %.15Lf\n", foil->x[n], foil->y[n]);

          fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[n],
                  foil->y[n], 0.0L, 2 * label1);
          fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[n],
                  foil->y[n], delz, 2 * label1 + 1);

          WriteToPoint(&foil->x[n], &foil->y[n], &(long double){0.0L},
                       &(int){2 * label1}, pointoutput[3]);
          WriteToPoint(&foil->x[n], &foil->y[n], &(long double){delz},
                       &(int){2 * label1 + 1}, pointoutput[3]);

          fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[n],
                  foil->y[n], 0.0L);
          fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[n],
                  foil->y[n], delz);
        };
      };
      if (foil->section[n] == 2) {
        if (foil->slope[n] > 0.0) {
          foil->x[n] -= cpar * cosl(atanl(foil->slope[n]));
          foil->y[n] -= cpar * sinl(atanl(foil->slope[n]));
        } else {
          foil->x[n] += cpar * cosl(atanl(foil->slope[n]));
          foil->y[n] += cpar * sinl(atanl(foil->slope[n]));
        };
        fprintf(pointoutput[0], "%.15Lf \t %.15Lf\n", foil->x[n], foil->y[n]);

        fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[n],
                foil->y[n], 0.0L, 2 * label1);
        fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[n],
                foil->y[n], delz, 2 * label1 + 1);

        WriteToPoint(&foil->x[n], &foil->y[n], &(long double){0.0L},
                     &(int){2 * label1}, pointoutput[3]);
        WriteToPoint(&foil->x[n], &foil->y[n], &(long double){delz},
                     &(int){2 * label1 + 1}, pointoutput[3]);

        fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[n],
                foil->y[n], 0.0L);
        fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[n],
                foil->y[n], delz);
      };

      n++;
      label1++;
    };
    n = 0;
    if (j <= 15) {
      cpar *= 1.03;

    } else {
      cpar += 0.08 * cpar;
    }

    j++;
  };
  int nodes;
  long double ynet[foil->count + 1];

  for (int i = 0; i < foil->count; i++) {
    if (i == 0) {
      printf("\nBefore integer division %.15Lf\n", (xL - foil->x[i]) / cpar);
      nodes = (xL - foil->x[i]) / cpar;
      printf("\nAfter integer division %d\n", (nodes));
      printf("\nAfter integer division %d\n", (nodes));
    };
    if (0 <= i && i <= ncrit + mp) {
      ynet[i] = foil->slope[i] * (xL - foil->x[i]) + foil->y[i];
    };
    if ((ncrit + mp < i) && (i < ntip)) {
      ynet[i] = (1.0 / foil->slope[i]) * (yU - foil->y[i]) + foil->x[i];
    };
    if (i == ntip) {
      ynet[i] = (1.0 / foil->slope[i]) * (yU - foil->y[i]) + foil->x[i];

      ynet[i + 1] = (-1.0 / foil->slope[i]) * (yB - ytail) + foil->x[i];
    };
    if ((i > ntip) && (i < foil->count - ncrit - (int)(r_num * mp / r_dem)) &&
        (foil->section[i] == 4)) {
      long double p1 = (foil->associatedT[i] - a.c * a.p) /
                       (foil->associatedT[ntip] - a.c * a.p);
      long double p2 = 1.0 - p1;
      long double tempslope =
          (0.0 - foil->slope[ntip]) * p1 + p2 * foil->slope[i];

      ynet[i + 1] = (1.0 / tempslope) * (yB - foil->y[i]) + foil->x[i];
    };
    if ((i > ntip) && (i < foil->count - ncrit - (int)(r_num * mp / r_dem)) &&
        (foil->section[i] != 4)) {
      ynet[i + 1] = (1.0 / foil->slope[i]) * (yB - foil->y[i]) + foil->x[i];
    };

    if (i >= foil->count - ncrit - (int)(r_num * mp / r_dem)) {

      ynet[i + 1] = foil->slope[i] * (xL - foil->x[i]) + foil->y[i];
    };
  };
  long double dx[foil->count + 1];
  long double dy[foil->count + 1];
  printf("\n NODES = %d \n", nodes);
  nodes *= -1;
  for (int u = 1; u <= nodes; u++) {

    for (int v = 0; v < foil->count; v++) {

      if (0 <= v && v <= ncrit + mp) {
        dx[v] = u * ((xL - foil->x[v]) / ((long double)(nodes)));
        dy[v] = u * ((ynet[v] - foil->y[v]) / ((long double)(nodes)));
        fprintf(pointoutput[0], "%.15Lf\t%.15Lf\n", foil->x[v] + dx[v],
                foil->y[v] + dy[v]);

        fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                foil->x[v] + dx[v], foil->y[v] + dy[v], 0.0L, 2 * label1);
        fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                foil->x[v] + dx[v], foil->y[v] + dy[v], delz, 2 * label1 + 1);

        WriteToPoint(&(long double){foil->x[v] + dx[v]},
                     &(long double){foil->y[v] + dy[v]}, &(long double){0.0L},
                     &(int){2 * label1}, pointoutput[3]);
        WriteToPoint(&(long double){foil->x[v] + dx[v]},
                     &(long double){foil->y[v] + dy[v]}, &(long double){delz},
                     &(int){2 * label1 + 1}, pointoutput[3]);

        fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[v] + dx[v],
                foil->y[v] + dy[v], 0.0L);
        fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[v] + dx[v],
                foil->y[v] + dy[v], delz);
      };
      if ((ncrit + mp < v) && (v < ntip)) {
        dx[v] = u * ((ynet[v] - foil->x[v]) / ((long double)(nodes)));
        dy[v] = u * ((yU - foil->y[v]) / ((long double)(nodes)));
        fprintf(pointoutput[0], "%.15Lf\t%.15Lf\n", foil->x[v] + dx[v],
                foil->y[v] + dy[v]);

        fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                foil->x[v] + dx[v], foil->y[v] + dy[v], 0.0L, 2 * label1);
        fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                foil->x[v] + dx[v], foil->y[v] + dy[v], delz, 2 * label1 + 1);

        WriteToPoint(&(long double){foil->x[v] + dx[v]},
                     &(long double){foil->y[v] + dy[v]}, &(long double){0.0L},
                     &(int){2 * label1}, pointoutput[3]);
        WriteToPoint(&(long double){foil->x[v] + dx[v]},
                     &(long double){foil->y[v] + dy[v]}, &(long double){delz},
                     &(int){2 * label1 + 1}, pointoutput[3]);

        fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[v] + dx[v],
                foil->y[v] + dy[v], 0.0L);
        fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[v] + dx[v],
                foil->y[v] + dy[v], delz);
      };
      if (v == ntip) {
        dx[v] = u * ((ynet[v] - foil->x[v]) / ((long double)(nodes)));
        dy[v] = u * ((yU - foil->y[v]) / ((long double)(nodes)));
        fprintf(pointoutput[0], "%.15Lf\t%.15Lf\n", foil->x[v] + dx[v],
                foil->y[v] + dy[v]);

        fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                foil->x[v] + dx[v], foil->y[v] + dy[v], 0.0L, 2 * label1);
        fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                foil->x[v] + dx[v], foil->y[v] + dy[v], delz, 2 * label1 + 1);

        WriteToPoint(&(long double){foil->x[v] + dx[v]},
                     &(long double){foil->y[v] + dy[v]}, &(long double){0.0L},
                     &(int){2 * label1}, pointoutput[3]);
        WriteToPoint(&(long double){foil->x[v] + dx[v]},
                     &(long double){foil->y[v] + dy[v]}, &(long double){delz},
                     &(int){2 * label1 + 1}, pointoutput[3]);

        fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[v] + dx[v],
                foil->y[v] + dy[v], 0.0L);
        fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[v] + dx[v],
                foil->y[v] + dy[v], delz);

        label1++;
        tracecount++;
        for (int i = 1; i <= backtrace; i++) {
          long double deltax =
              (xR - (foil->x[v] + dx[v])) / ((long double)(backtrace));
          trace1->x[i - 1] = ((long double)(i)) * deltax + foil->x[v] + dx[v];
          trace1->y = foil->y[v] + dy[v];
          trace1->label[i - 1] = label1;
          fprintf(pointoutput[0], "%.15Lf\t%.15Lf\n", trace1->x[i - 1],
                  trace1->y);

          fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                  trace1->x[i - 1], trace1->y, 0.0L, 2 * label1);
          fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                  trace1->x[i - 1], trace1->y, delz, 2 * label1 + 1);

          WriteToPoint(&trace1->x[i - 1], &trace1->y, &(long double){0.0L},
                       &(int){2 * label1}, pointoutput[3]);
          WriteToPoint(&trace1->x[i - 1], &trace1->y, &(long double){delz},
                       &(int){2 * label1 + 1}, pointoutput[3]);

          fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", trace1->x[i - 1],
                  trace1->y, 0.0L);
          fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", trace1->x[i - 1],
                  trace1->y, delz);

          label1++;
        };

        dx[v + 1] = u * ((ynet[v + 1] - foil->x[v]) / ((long double)(nodes)));
        dy[v + 1] = u * ((yB - ytail) / ((long double)(nodes)));
        fprintf(pointoutput[0], "%.15Lf\t%.15Lf\n", foil->x[v] + dx[v + 1],
                ytail + dy[v + 1]);

        fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                foil->x[v] + dx[v + 1], ytail + dy[v + 1], 0.0L, 2 * label1);
        fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                foil->x[v] + dx[v + 1], ytail + dy[v + 1], delz,
                2 * label1 + 1);

        WriteToPoint(&(long double){foil->x[v] + dx[v + 1]},
                     &(long double){ytail + dy[v + 1]}, &(long double){0.0L},
                     &(int){2 * label1}, pointoutput[3]);
        WriteToPoint(&(long double){foil->x[v] + dx[v + 1]},
                     &(long double){ytail + dy[v + 1]}, &(long double){delz},
                     &(int){2 * label1 + 1}, pointoutput[3]);

        fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n",
                foil->x[v] + dx[v + 1], ytail + dy[v + 1], 0.0L);
        fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n",
                foil->x[v] + dx[v + 1], ytail + dy[v + 1], delz);

        tracecount++;
        label1++;
        for (int i = 1; i <= backtrace; i++) {
          long double deltax =
              (xR - (foil->x[v] + dx[v + 1])) / ((long double)(backtrace));
          trace1->x[i - 1] =
              ((long double)(i)) * deltax + foil->x[v] + dx[v + 1];
          trace1->y = ytail + dy[v + 1];
          trace1->label[i - 1] = label1;
          fprintf(pointoutput[0], "%.15Lf\t%.15Lf\n", trace1->x[i - 1],
                  trace1->y);

          fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                  trace1->x[i - 1], trace1->y, 0.0L, 2 * label1);
          fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                  trace1->x[i - 1], trace1->y, delz, 2 * label1 + 1);

          WriteToPoint(&trace1->x[i - 1], &trace1->y, &(long double){0.0L},
                       &(int){2 * label1}, pointoutput[3]);
          WriteToPoint(&trace1->x[i - 1], &trace1->y, &(long double){delz},
                       &(int){2 * label1 + 1}, pointoutput[3]);

          fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", trace1->x[i - 1],
                  trace1->y, 0.0L);
          fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", trace1->x[i - 1],
                  trace1->y, delz);

          label1++;
        };
        label1--;
      };
      if ((v > ntip) && (v < foil->count - ncrit - (int)(r_num * mp / r_dem)) &&
          (foil->section[v] == 4)) {
        dx[v + 1] = u * ((ynet[v + 1] - foil->x[v]) / ((long double)(nodes)));
        dy[v + 1] = u * ((yB - foil->y[v]) / ((long double)(nodes)));

        fprintf(pointoutput[0], "%.15Lf\t%.15Lf\n", foil->x[v] + dx[v + 1],
                foil->y[v] + dy[v + 1]);

        fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                foil->x[v] + dx[v + 1], foil->y[v] + dy[v + 1], 0.0L,
                2 * label1);
        fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                foil->x[v] + dx[v + 1], foil->y[v] + dy[v + 1], delz,
                2 * label1 + 1);

        WriteToPoint(&(long double){foil->x[v] + dx[v + 1]},
                     &(long double){foil->y[v] + dy[v + 1]},
                     &(long double){0.0L}, &(int){2 * label1}, pointoutput[3]);
        WriteToPoint(&(long double){foil->x[v] + dx[v + 1]},
                     &(long double){foil->y[v] + dy[v + 1]},
                     &(long double){delz}, &(int){2 * label1 + 1},
                     pointoutput[3]);

        fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n",
                foil->x[v] + dx[v + 1], foil->y[v] + dy[v + 1], 0.0L);
        fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n",
                foil->x[v] + dx[v + 1], foil->y[v] + dy[v + 1], delz);
      }

      if ((v > ntip) && (v < foil->count - ncrit - (int)(r_num * mp / r_dem)) &&
          (foil->section[v] != 4)) {
        dx[v + 1] = u * ((ynet[v + 1] - foil->x[v]) / ((long double)(nodes)));
        dy[v + 1] = u * ((yB - foil->y[v]) / ((long double)(nodes)));
        fprintf(pointoutput[0], "%.15Lf\t%.15Lf\n", foil->x[v] + dx[v + 1],
                foil->y[v] + dy[v + 1]);

        fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                foil->x[v] + dx[v + 1], foil->y[v] + dy[v + 1], 0.0L,
                2 * label1);
        fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                foil->x[v] + dx[v + 1], foil->y[v] + dy[v + 1], delz,
                2 * label1 + 1);

        WriteToPoint(&(long double){foil->x[v] + dx[v + 1]},
                     &(long double){foil->y[v] + dy[v + 1]},
                     &(long double){0.0L}, &(int){2 * label1}, pointoutput[3]);
        WriteToPoint(&(long double){foil->x[v] + dx[v + 1]},
                     &(long double){foil->y[v] + dy[v + 1]},
                     &(long double){delz}, &(int){2 * label1 + 1},
                     pointoutput[3]);

        fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n",
                foil->x[v] + dx[v + 1], foil->y[v] + dy[v + 1], 0.0L);
        fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n",
                foil->x[v] + dx[v + 1], foil->y[v] + dy[v + 1], delz);
      };
      if (v >= foil->count - ncrit - (int)(r_num * mp / r_dem)) {
        dx[v + 1] = u * ((xL - foil->x[v]) / ((long double)(nodes)));
        dy[v + 1] = u * ((ynet[v + 1] - foil->y[v]) / ((long double)(nodes)));

        fprintf(pointoutput[0], "%.15Lf\t%.15Lf\n", foil->x[v] + dx[v + 1],
                foil->y[v] + dy[v + 1]);

        fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                foil->x[v] + dx[v + 1], foil->y[v] + dy[v + 1], 0.0L,
                2 * label1);
        fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n",
                foil->x[v] + dx[v + 1], foil->y[v] + dy[v + 1], delz,
                2 * label1 + 1);

        WriteToPoint(&(long double){foil->x[v] + dx[v + 1]},
                     &(long double){foil->y[v] + dy[v + 1]},
                     &(long double){0.0L}, &(int){2 * label1}, pointoutput[3]);
        WriteToPoint(&(long double){foil->x[v] + dx[v + 1]},
                     &(long double){foil->y[v] + dy[v + 1]},
                     &(long double){delz}, &(int){2 * label1 + 1},
                     pointoutput[3]);

        fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n",
                foil->x[v] + dx[v + 1], foil->y[v] + dy[v + 1], 0.0L);
        fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n",
                foil->x[v] + dx[v + 1], foil->y[v] + dy[v + 1], delz);
      };
      label1++;
    };
  };
  // fgets(buffer, sizeof(buffer), dafile);
  fprintf(pointoutput[2], "\n )\n");
  fclose(pointoutput[0]);
  fclose(pointoutput[1]);
  fclose(pointoutput[2]);
  fclose(pointoutput[3]);
  // pointoutput[2] = fopen("points", "r");
  faces[0] = fopen("faces", "a");
  faces[1] = fopen("fitago", "w");
  faces[2] = fopen("facecheck", "w");
  owner[0] = fopen("owner", "a");
  owner[1] = fopen("owncheck", "w");

  neighbour[0] = fopen("neighbour", "a");
  neighbour[1] = fopen("neighbourcheck", "w");

  int bc = backtrace;

  int nFaces[8];

  nFaces[FrontAndBack] = 2 * (2 * ntip + 2 * bc) * (j + nodes - 1);
  nFaces[Aerofoil] = 2 * ntip;
  nFaces[Inlet] = mp + (mp * r_num) / r_dem + 2 * ncrit;
  nFaces[Outlet] = 2 * (j + nodes - 1);
  nFaces[Top1] = ntip - mp - ncrit;
  nFaces[Bottom1] = ntip - ncrit - (mp * r_num) / r_dem;
  nFaces[Top2] = bc;
  nFaces[Bottom2] = bc;

  int n_InternalFaces = 2 * (2 * ntip + 2 * bc) * (j + nodes - 1) +
                        (2 * ntip + bc) +
                        (2 * ntip + 2 * bc) * (j + nodes - 1) +
                        (2 * ntip + 1 + 2 * bc) * (j + nodes - 1) -
                        (nFaces[0] + nFaces[1] + nFaces[2] + nFaces[3] +
                         nFaces[4] + nFaces[5] + nFaces[6] + nFaces[7]);

  int startFace[8];
  // startFace[FrontAndBack] = n_InternalFaces;
  startFace[Aerofoil] = n_InternalFaces + nFaces[FrontAndBack];
  startFace[Inlet] = n_InternalFaces + nFaces[FrontAndBack] + nFaces[Aerofoil];
  startFace[Outlet] = startFace[Inlet] + nFaces[Inlet];
  startFace[Top1] = startFace[Outlet] + nFaces[Outlet];
  startFace[Bottom1] = startFace[Top1] + nFaces[Top1];
  startFace[Top2] = startFace[Bottom1] + nFaces[Bottom1];
  startFace[Bottom2] = startFace[Top2] + nFaces[Top2];

  fprintf(faces[0], "\n\n");
  fprintf(faces[0], "%d",
          2 * (2 * ntip + 2 * bc) * (j + nodes - 1) + (2 * ntip + bc) +
              (2 * ntip + 2 * bc) * (j + nodes - 1) +
              (2 * ntip + 1 + 2 * bc) * (j + nodes - 1));
  fprintf(faces[0], "\n(\n");
  fprintf(owner[0], "\n\n%d\n(\n",
          2 * (2 * ntip + 2 * bc) * (j + nodes - 1) + (2 * ntip + bc) +
              (2 * ntip + 2 * bc) * (j + nodes - 1) +
              (2 * ntip + 1 + 2 * bc) * (j + nodes - 1));

  fprintf(neighbour[0], "\n\n%d\n(\n",
          2 * (2 * ntip + 2 * bc) * (j + nodes - 1) + (2 * ntip + bc) +
              (2 * ntip + 2 * bc) * (j + nodes - 1) +
              (2 * ntip + 1 + 2 * bc) * (j + nodes - 1) -
              (nFaces[0] + nFaces[1] + nFaces[2] + nFaces[3] + nFaces[4] +
               nFaces[5] + nFaces[6] + nFaces[7]));

  int loop1 = 2 * (2 * ntip + backtrace);
  int loop2 = 2 * (2 * ntip + 2 * backtrace + 1);

  int flabel = 0;
  // typedef enum { AEROFOIL, NORTH, SOUTH, EAST, WEST, IN, OUT, INLET, OUTLET,
  // TOP, BOTTOM, REVERSE, TOP2, BOTTOM2 } facecode;

  facecode bcode;

  int t2n =
      2 * ntip + 2 * backtrace + 2 * ntip - 2 * ncrit -
      2 * ((int)((((long double)(r_num)) / ((long double)(r_dem))) * (mp)));
  bcode = AEROFOIL;
  int total = loop1 + (j + nodes - 1) * loop2;
  int sik = (j + nodes - 1) * (2 * ntip + 2 * backtrace);

  cell *oker;
  oker = (cell *)malloc(sik * sizeof(cell));

  int i1, comp;
  int k = 0;
  int clabel = 0;

  int AerofoilF[nFaces[Aerofoil]];
  int InletF[nFaces[Inlet]];
  int OutletF[nFaces[Outlet]];
  int Top1F[nFaces[Top1]];
  int Bottom1F[nFaces[Bottom1]];
  int Top2F[nFaces[Top2]];
  int Bottom2F[nFaces[Bottom2]];

  // Assigning Boundary Face numbers

  // Aerofoil
  bcode = AEROFOIL;
  flabel = startFace[Aerofoil];

  for (int i = 0; i < nFaces[Aerofoil]; i++) {
    AerofoilF[i] = flabel;
    flabel++;
  };

  // Inlet

  bcode = INLET;

  for (int i = 0; i < nFaces[Inlet]; i++) {
    InletF[i] = flabel;
    flabel++;
  };
  // Printing Outlet
  bcode = OUTLET;

  for (int i = 0; i < nFaces[Outlet]; i++) {
    OutletF[i] = flabel;
    flabel++;
  };
  // Top1
  bcode = TOP;

  for (int i = 0; i < nFaces[Top1]; i++) {
    Top1F[i] = flabel;
    flabel++;
  };

  // Bottom1
  bcode = BOTTOM;

  for (int i = 0; i < nFaces[Bottom1]; i++) {
    Bottom1F[i] = flabel;
    flabel++;
  };
  // Top2
  bcode = TOP2;

  for (int i = 0; i < nFaces[Top2]; i++) {
    Top2F[i] = flabel;
    flabel++;
  };
  // Bottom2
  bcode = BOTTOM2;

  for (int i = 0; i < nFaces[Bottom2]; i++) {
    Bottom2F[i] = flabel;
    flabel++;
  };

  flabel = 0;
  // TailTip Printing
  clabel = 2 * ntip;
  bcode = NORTH;
  int Tippy[bc];

  printf("\nTippy Declared succesfully with size %d\n", bc);

  int P55 = 0;

  for (int i = 2 * ntip; i < 2 * ntip + 2 * bc; i += 2) {
    ownFILL(owner, neighbour, clabel, clabel + bc, &bcode);
    onlyFILL(faces, i, i + 2, i + 3, i + 1, &flabel, bcode);
    Tippy[P55] = flabel;
    flabel++;
    clabel++;
    P55++;
  }

  printf("\nTippy correct\n");
  clabel = 0;

  bcode = NORTH;

  int I = 0;

  for (int i = t2n + loop2; i < loop1 + loop2 - 2; i += 2) {
    oker[clabel].faces[0] = AerofoilF[I];
    ownFILL(owner, neighbour, clabel, clabel + 2 * ntip + 2 * bc, &bcode);
    faceFILL(faces, i, i + 1, i + 3, i + 2, &flabel, bcode, &oker[clabel], 1,
             &clabel);
    I++;
  };
  i1 = loop1 + loop2 - 2;
  comp = loop1;
  oker[clabel].faces[0] = AerofoilF[I];
  ownFILL(owner, neighbour, clabel, clabel + 2 * ntip + 2 * bc, &bcode);
  faceFILL(faces, i1, i1 + 1, comp + 1, comp, &flabel, bcode, &oker[clabel], 1,
           &clabel);
  I++;

  for (int i = loop1; i < loop1 + 2 * ntip; i += 2) {
    oker[clabel].faces[0] = AerofoilF[I];
    ownFILL(owner, neighbour, clabel, clabel + 2 * ntip + 2 * bc, &bcode);
    faceFILL(faces, i, i + 1, i + 3, i + 2, &flabel, bcode, &oker[clabel], 1,
             &clabel);
    I++;
  };
  i1 = 2 * ntip + loop2;
  comp = i1 + 2 * bc + 2;
  oker[clabel].faces[0] = AerofoilF[I];
  ownFILL(owner, neighbour, clabel, clabel + 2 * ntip + 2 * bc, &bcode);
  faceFILL(faces, i1, i1 + 1, comp + 1, comp, &flabel, bcode, &oker[clabel], 1,
           &clabel);
  I++;

  for (int i = comp; i < t2n + loop2; i += 2) {
    oker[clabel].faces[0] = AerofoilF[I];
    ownFILL(owner, neighbour, clabel, clabel + 2 * ntip + 2 * bc, &bcode);
    faceFILL(faces, i, i + 1, i + 3, i + 2, &flabel, bcode, &oker[clabel], 1,
             &clabel);
    I++;
  };

  I = 0;

  for (int i = 2 * ntip + loop1; i < 2 * ntip + loop1 + 2 * bc; i += 2) {
    oker[clabel].faces[0] = Tippy[I];
    ownFILL(owner, neighbour, clabel, clabel + 2 * ntip + 2 * bc, &bcode);
    faceFILL(faces, i, i + 1, i + 3, i + 2, &flabel, bcode, &oker[clabel], 1,
             &clabel);
    I++;
  }

  I = 0;

  for (int i = 2 * ntip + loop2; i < 2 * ntip + loop2 + 2 * bc; i += 2) {
    oker[clabel].faces[0] = Tippy[I];
    ownFILL(owner, neighbour, clabel, clabel + 2 * ntip + 2 * bc, &bcode);
    faceFILL(faces, i, i + 2, i + 3, i + 1, &flabel, bcode, &oker[clabel], 1,
             &clabel);
    I++;
  };

  k = 2;

  while (k < j + nodes - 1) {

    for (int i = t2n + loop2 * k; i < loop1 + k * loop2 - 2; i += 2) {
      oker[clabel].faces[0] = flabel - 2 * ntip - 2 * bc;
      ownFILL(owner, neighbour, clabel, clabel + 2 * ntip + 2 * bc, &bcode);
      faceFILL(faces, i, i + 1, i + 3, i + 2, &flabel, bcode, &oker[clabel], 1,
               &clabel);
    };

    i1 = loop1 + k * loop2 - 2;
    comp = loop1 + (k - 1) * loop2;
    oker[clabel].faces[0] = flabel - 2 * ntip - 2 * bc;
    ownFILL(owner, neighbour, clabel, clabel + 2 * ntip + 2 * bc, &bcode);
    faceFILL(faces, i1, i1 + 1, comp + 1, comp, &flabel, bcode, &oker[clabel],
             1, &clabel);

    for (int i = loop1 + (k - 1) * loop2;
         i < loop1 + (k - 1) * loop2 + 2 * ntip; i += 2) {
      oker[clabel].faces[0] = flabel - 2 * ntip - 2 * bc;
      ownFILL(owner, neighbour, clabel, clabel + 2 * ntip + 2 * bc, &bcode);
      faceFILL(faces, i, i + 1, i + 3, i + 2, &flabel, bcode, &oker[clabel], 1,
               &clabel);
    }

    i1 = 2 * ntip + k * loop2;
    comp = i1 + 2 * bc + 2;
    oker[clabel].faces[0] = flabel - 2 * ntip - 2 * bc;
    ownFILL(owner, neighbour, clabel, clabel + 2 * ntip + 2 * bc, &bcode);
    faceFILL(faces, i1, i1 + 1, comp + 1, comp, &flabel, bcode, &oker[clabel],
             1, &clabel);

    for (int i = comp; i < t2n + loop2 * k; i += 2) {
      oker[clabel].faces[0] = flabel - 2 * ntip - 2 * bc;
      ownFILL(owner, neighbour, clabel, clabel + 2 * ntip + 2 * bc, &bcode);
      faceFILL(faces, i, i + 1, i + 3, i + 2, &flabel, bcode, &oker[clabel], 1,
               &clabel);
    };

    for (int i = loop1 + (k - 1) * loop2 + 2 * ntip;
         i < loop1 + (k - 1) * loop2 + 2 * ntip + 2 * bc; i += 2) {
      oker[clabel].faces[0] = flabel - 2 * ntip - 2 * bc;
      ownFILL(owner, neighbour, clabel, clabel + 2 * ntip + 2 * bc, &bcode);
      faceFILL(faces, i, i + 1, i + 3, i + 2, &flabel, bcode, &oker[clabel], 1,
               &clabel);
    };

    for (int i = 2 * ntip + k * loop2; i < 2 * ntip + k * loop2 + 2 * bc;
         i += 2) {
      oker[clabel].faces[0] = flabel - 2 * ntip - 2 * bc;
      ownFILL(owner, neighbour, clabel, clabel + 2 * ntip + 2 * bc, &bcode);
      faceFILL(faces, i, i + 2, i + 3, i + 1, &flabel, bcode, &oker[clabel], 1,
               &clabel);
    };

    k++;
  };

  int dlabel = flabel - 2 * ntip - 2 * bc;

  for (int i = 0; i < nFaces[Inlet]; i++) {
    oker[clabel].faces[0] = dlabel;
    oker[clabel].faces[1] = InletF[i];
    dlabel++;
    clabel++;
  };

  for (int i = 0; i < nFaces[Top1]; i++) {
    oker[clabel].faces[0] = dlabel;
    oker[clabel].faces[1] = Top1F[i];
    dlabel++;
    clabel++;
  };

  for (int i = 0; i < nFaces[Bottom1]; i++) {
    oker[clabel].faces[0] = dlabel;
    oker[clabel].faces[1] = Bottom1F[i];
    dlabel++;
    clabel++;
  };

  for (int i = 0; i < nFaces[Top2]; i++) {
    oker[clabel].faces[0] = dlabel;
    oker[clabel].faces[1] = Top2F[i];
    dlabel++;
    clabel++;
  };

  for (int i = 0; i < nFaces[Bottom2]; i++) {
    oker[clabel].faces[0] = dlabel;
    oker[clabel].faces[1] = Bottom2F[i];
    dlabel++;
    clabel++;
  };

  printf("North faces printed correctly\n");

  // Printing Verticial Faces
  bcode = WEST;

  clabel = 0;
  i1 = t2n;
  comp = t2n + loop2;
  I = flabel;
  ownFILL(owner, neighbour, clabel, clabel + 2 * ntip - 1, &bcode);
  oker[clabel].faces[2] = flabel + 1;
  faceFILL(faces, i1, i1 + 1, comp + 1, comp, &flabel, bcode, &oker[clabel], 3,
           &clabel);

  for (int i = t2n + 2; i < loop1; i += 2) {
    ownFILL(owner, neighbour, clabel - 1, clabel, &bcode);
    oker[clabel].faces[2] = flabel + 1;
    faceFILL(faces, i, i + loop2, i + loop2 + 1, i + 1, &flabel, bcode,
             &oker[clabel], 3, &clabel);
  };

  for (int i = 0; i < 2 * ntip; i += 2) {
    ownFILL(owner, neighbour, clabel - 1, clabel, &bcode);
    oker[clabel].faces[2] = flabel + 1;
    faceFILL(faces, i, i + loop1, i + loop1 + 1, i + 1, &flabel, bcode,
             &oker[clabel], 3, &clabel);
  };
  dlabel = flabel;
  ownFILL(owner, neighbour, clabel - 1, clabel - 1 + nFaces[Bottom1] + 1,
          &bcode);
  onlyFILL(faces, 2 * ntip, 2 * ntip + loop1, 2 * ntip + loop1 + 1,
           2 * ntip + 1, &flabel, bcode);
  flabel++;

  oker[clabel].faces[2] = flabel + 1;
  ownFILL(owner, neighbour, clabel, clabel + nFaces[Bottom1] + bc, &bcode);
  faceFILL(faces, 2 * ntip, 2 * ntip + 1, 2 * ntip + 1 + loop2,
           2 * ntip + loop2, &flabel, bcode, &oker[clabel], 3, &clabel);

  for (int i = 2 * ntip + 2 * bc + 2; i < t2n - 2; i += 2) {
    oker[clabel].faces[2] = flabel + 1;
    ownFILL(owner, neighbour, clabel - 1, clabel, &bcode);
    faceFILL(faces, i, i + loop2, i + loop2 + 1, i + 1, &flabel, bcode,
             &oker[clabel], 3, &clabel);
  };

  oker[clabel].faces[2] = I;
  ownFILL(owner, neighbour, clabel - 1, clabel, &bcode);
  i1 = t2n - 2;
  faceFILL(faces, i1, i1 + loop2, i1 + loop2 + 1, i1 + 1, &flabel, bcode,
           &oker[clabel], 3, &clabel);

  i1 = 2 * ntip + 2;
  oker[clabel].faces[3] = dlabel;
  ownFILL(owner, neighbour, clabel, clabel + 1, &bcode);
  faceFILL(faces, i1, i1 + loop1, i1 + loop1 + 1, i1 + 1, &flabel, bcode,
           &oker[clabel], 2, &clabel);

  for (int i = 2 * ntip + 4; i < 2 * ntip + 2 * bc; i += 2) {
    oker[clabel].faces[3] = flabel - 1;
    ownFILL(owner, neighbour, clabel, clabel + 1, &bcode);
    faceFILL(faces, i, i + loop1, i + loop1 + 1, i + 1, &flabel, bcode,
             &oker[clabel], 2, &clabel);
  };
  oker[clabel].faces[3] = flabel - 1;
  oker[clabel].faces[2] = OutletF[j + nodes - 1];
  clabel++;

  oker[clabel].faces[3] = dlabel + 1;
  ownFILL(owner, neighbour, clabel, clabel + 1, &bcode);
  faceFILL(faces, i1, i1 + 1, i1 + 1 + loop2, i1 + loop2, &flabel, bcode,
           &oker[clabel], 2, &clabel);

  for (int i = 2 * ntip + 4; i < 2 * ntip + 2 * bc; i += 2) {
    oker[clabel].faces[3] = flabel - 1;
    ownFILL(owner, neighbour, clabel, clabel + 1, &bcode);
    faceFILL(faces, i, i + 1, i + 1 + loop2, i + loop2, &flabel, bcode,
             &oker[clabel], 2, &clabel);
  };
  oker[clabel].faces[3] = flabel - 1;
  oker[clabel].faces[2] = OutletF[j + nodes - 2];
  clabel++;

  printf("\nclabel is supposed to be %d, but its %d\n", 2 * ntip + 2 * bc,
         clabel);
  k = 1;

  while (k < j + nodes - 1) {

    i1 = t2n + k * loop2;
    comp = i1 + loop2;
    I = flabel;
    ownFILL(owner, neighbour, clabel, clabel + 2 * ntip - 1, &bcode);
    oker[clabel].faces[2] = flabel + 1;
    faceFILL(faces, i1, i1 + 1, comp + 1, comp, &flabel, bcode, &oker[clabel],
             3, &clabel);
    for (int i = t2n + k * loop2 + 2; i < loop1 + k * loop2; i += 2) {
      ownFILL(owner, neighbour, clabel - 1, clabel, &bcode);
      oker[clabel].faces[2] = flabel + 1;
      faceFILL(faces, i, i + loop2, i + loop2 + 1, i + 1, &flabel, bcode,
               &oker[clabel], 3, &clabel);
    }

    for (int i = loop1 + (k - 1) * loop2;
         i < loop1 + (k - 1) * loop2 + 2 * ntip; i += 2) {
      ownFILL(owner, neighbour, clabel - 1, clabel, &bcode);
      //     printf("\nOker size = %d, printing to %d", sik, clabel);
      oker[clabel].faces[2] = flabel + 1;
      faceFILL(faces, i, i + loop2, i + loop2 + 1, i + 1, &flabel, bcode,
               &oker[clabel], 3, &clabel);
    }

    dlabel = flabel;
    i1 = loop1 + (k - 1) * loop2 + 2 * ntip;
    comp = i1 + loop2;
    ownFILL(owner, neighbour, clabel - 1, clabel - 1 + nFaces[Bottom1] + 1,
            &bcode);
    onlyFILL(faces, i1, comp, comp + 1, i1 + 1, &flabel, bcode);
    flabel++;

    oker[clabel].faces[2] = flabel + 1;
    ownFILL(owner, neighbour, clabel, clabel + nFaces[Bottom1] + bc, &bcode);
    i1 = 2 * ntip + k * loop2;
    comp = i1 + loop2;
    faceFILL(faces, i1, i1 + 1, comp + 1, comp, &flabel, bcode, &oker[clabel],
             3, &clabel);

    for (int i = 2 * ntip + k * loop2 + 2 * bc + 2; i < t2n + k * loop2 - 2;
         i += 2) {
      ownFILL(owner, neighbour, clabel - 1, clabel, &bcode);
      oker[clabel].faces[2] = flabel + 1;
      faceFILL(faces, i, i + loop2, i + loop2 + 1, i + 1, &flabel, bcode,
               &oker[clabel], 3, &clabel);
    };

    i1 = t2n + k * loop2 - 2;
    comp = i1 + loop2;
    ownFILL(owner, neighbour, clabel - 1, clabel, &bcode);
    oker[clabel].faces[2] = I;
    faceFILL(faces, i1, comp, comp + 1, i1 + 1, &flabel, bcode, &oker[clabel],
             3, &clabel);

    i1 = 2 * ntip + 2 + loop1 + (k - 1) * loop2;
    comp = i1 + loop2;
    ownFILL(owner, neighbour, clabel, clabel + 1, &bcode);
    oker[clabel].faces[3] = dlabel;
    faceFILL(faces, i1, comp, comp + 1, i1 + 1, &flabel, bcode, &oker[clabel],
             2, &clabel);

    for (int i = 2 * ntip + loop1 + (k - 1) * loop2 + 4;
         i < 2 * ntip + loop1 + (k - 1) * loop2 + 2 * bc; i += 2) {

      ownFILL(owner, neighbour, clabel, clabel + 1, &bcode);
      oker[clabel].faces[3] = flabel - 1;
      faceFILL(faces, i, i + loop2, i + loop2 + 1, i + 1, &flabel, bcode,
               &oker[clabel], 2, &clabel);
    };
    oker[clabel].faces[3] = flabel - 1;
    oker[clabel].faces[2] = OutletF[j + nodes - 1 + k];
    clabel++;

    oker[clabel].faces[3] = dlabel + 1;
    ownFILL(owner, neighbour, clabel, clabel + 1, &bcode);
    i1 = 2 * ntip + 2 + k * loop2;
    comp = i1 + loop2;
    faceFILL(faces, i1, i1 + 1, comp + 1, comp, &flabel, bcode, &oker[clabel],
             2, &clabel);

    for (int i = 2 * ntip + 4 + k * loop2; i < 2 * ntip + k * loop2 + 2 * bc;
         i += 2) {

      oker[clabel].faces[3] = flabel - 1;
      ownFILL(owner, neighbour, clabel, clabel + 1, &bcode);
      faceFILL(faces, i, i + 1, i + 1 + loop2, i + loop2, &flabel, bcode,
               &oker[clabel], 2, &clabel);
    }

    oker[clabel].faces[3] = flabel - 1;
    oker[clabel].faces[2] = OutletF[j + nodes - 2 - k];
    clabel++;

    k++;
  };

  printf("Vertical Faces printed correctly\n");

  // Boundary faces begin printing here
  clabel = 0;
  bcode = IN;

  startFace[FrontAndBack] = flabel + 0;

  for (int i = t2n; i < loop1 - 2; i++) {
    ownFILL(owner, neighbour, clabel, -1, &bcode);
    faceFILL(faces, i, i + loop2, i + loop2 + 2, i + 2, &flabel, bcode,
             &oker[clabel], 4, &clabel);

    i++;
  };
  i1 = loop1 - 2;
  comp = 0;
  ownFILL(owner, neighbour, clabel, -1, &bcode);
  faceFILL(faces, i1, i1 + loop2, loop1, 0, &flabel, bcode, &oker[clabel], 4,
           &clabel);

  for (int i = 0; i < 2 * ntip; i++) {
    ownFILL(owner, neighbour, clabel, -1, &bcode);
    faceFILL(faces, i, i + loop1, i + loop1 + 2, i + 2, &flabel, bcode,
             &oker[clabel], 4, &clabel);

    i++;
  };
  i1 = 2 * ntip;
  comp = 2 * ntip + 2 * bc + 2;
  ownFILL(owner, neighbour, clabel, -1, &bcode);
  faceFILL(faces, i1, i1 + loop2, comp + loop2, comp, &flabel, bcode,
           &oker[clabel], 4, &clabel);

  for (int i = 2 * ntip + 2 * bc + 2; i < t2n; i++) {
    ownFILL(owner, neighbour, clabel, -1, &bcode);
    faceFILL(faces, i, i + loop2, i + loop2 + 2, i + 2, &flabel, bcode,
             &oker[clabel], 4, &clabel);
    i++;
  };

  for (int i = 2 * ntip + 2; i <= 2 * ntip + 2 * bc; i++) {

    ownFILL(owner, neighbour, clabel, -1, &bcode);
    faceFILL(faces, i, i - 2, i - 2 + loop1, i + loop1, &flabel, bcode,
             &oker[clabel], 4, &clabel);
    i++;
  };

  for (int i = 2 * ntip + 2; i <= 2 * ntip + 2 * bc; i++) {

    ownFILL(owner, neighbour, clabel, -1, &bcode);
    faceFILL(faces, i, i + loop2, i + loop2 - 2, i - 2, &flabel, bcode,
             &oker[clabel], 4, &clabel);
    i++;
  };

  k = 1;
  while (k < j + nodes - 1) {

    for (int i = t2n + k * loop2; i < loop1 + k * loop2 - 2; i++) {
      ownFILL(owner, neighbour, clabel, -1, &bcode);
      faceFILL(faces, i, i + loop2, i + loop2 + 2, i + 2, &flabel, bcode,
               &oker[clabel], 4, &clabel);
      i++;
    };
    i1 = loop1 + k * loop2 - 2;
    comp = loop1 + (k - 1) * loop2;
    ownFILL(owner, neighbour, clabel, -1, &bcode);
    faceFILL(faces, i1, i1 + loop2, comp + loop2, comp, &flabel, bcode,
             &oker[clabel], 4, &clabel);

    for (int i = loop1 + (k - 1) * loop2;
         i < loop1 + (k - 1) * loop2 + 2 * ntip; i++) {

      ownFILL(owner, neighbour, clabel, -1, &bcode);
      faceFILL(faces, i, i + loop2, i + loop2 + 2, i + 2, &flabel, bcode,
               &oker[clabel], 4, &clabel);
      i++;
    };

    i1 = 2 * ntip + loop2 * k;
    comp = 2 * ntip + 2 * bc + 2 + loop2 * k;
    ownFILL(owner, neighbour, clabel, -1, &bcode);
    faceFILL(faces, i1, i1 + loop2, comp + loop2, comp, &flabel, bcode,
             &oker[clabel], 4, &clabel);

    for (int i = 2 * ntip + 2 * bc + 2 + loop2 * k; i < t2n + loop2 * k; i++) {

      ownFILL(owner, neighbour, clabel, -1, &bcode);
      faceFILL(faces, i, i + loop2, i + loop2 + 2, i + 2, &flabel, bcode,
               &oker[clabel], 4, &clabel);
      i++;
    };
    for (int i = 2 * ntip + 2 + loop1 + (k - 1) * loop2;
         i <= 2 * ntip + loop1 + (k - 1) * loop2 + 2 * bc; i++) {

      ownFILL(owner, neighbour, clabel, -1, &bcode);
      faceFILL(faces, i, i - 2, i - 2 + loop2, i + loop2, &flabel, bcode,
               &oker[clabel], 4, &clabel);
      i++;
    };

    for (int i = 2 * ntip + 4 + loop1 + (k - 1) * loop2 + 2 * bc;
         i <= 2 * ntip + 2 + loop1 + (k - 1) * loop2 + 4 * bc; i++) {

      ownFILL(owner, neighbour, clabel, -1, &bcode);
      faceFILL(faces, i, i + loop2, i + loop2 - 2, i - 2, &flabel, bcode,
               &oker[clabel], 4, &clabel);
      i++;
    };

    k++;
  };

  printf("\nLast cell printed to is (3) %d\n", clabel - 1);
  k = 1;
  clabel = 0;

  bcode = OUT;
  for (int i = t2n; i < loop1 - 2; i++) {
    ownFILL(owner, neighbour, clabel, -1, &bcode);
    faceFILL(faces, i + 3, i + loop2 + 3, i + loop2 + 1, i + 1, &flabel, bcode,
             &oker[clabel], 5, &clabel);
    i++;
  };
  i1 = loop1 - 2;
  comp = 0;

  ownFILL(owner, neighbour, clabel, -1, &bcode);
  GfaceFILL(faces, i1, i1 + loop2, loop1, 0, &flabel, bcode, &oker[clabel], 4,
            &clabel);

  for (int i = 0; i < 2 * ntip; i++) {
    ownFILL(owner, neighbour, clabel, -1, &bcode);
    GfaceFILL(faces, i, i + loop1, i + loop1 + 2, i + 2, &flabel, bcode,
              &oker[clabel], 4, &clabel);
    i++;
  };
  i1 = 2 * ntip;
  comp = 2 * ntip + 2 * bc + 2;

  ownFILL(owner, neighbour, clabel, -1, &bcode);
  GfaceFILL(faces, i1, i1 + loop2, comp + loop2, comp, &flabel, bcode,
            &oker[clabel], 4, &clabel);

  for (int i = 2 * ntip + 2 * bc + 2; i < t2n; i++) {
    ownFILL(owner, neighbour, clabel, -1, &bcode);
    GfaceFILL(faces, i, i + loop2, i + loop2 + 2, i + 2, &flabel, bcode,
              &oker[clabel], 4, &clabel);
    i++;
  };

  for (int i = 2 * ntip + 2; i <= 2 * ntip + 2 * bc; i++) {

    ownFILL(owner, neighbour, clabel, -1, &bcode);
    GfaceFILL(faces, i, i - 2, i - 2 + loop1, i + loop1, &flabel, bcode,
              &oker[clabel], 4, &clabel);
    i++;
  };

  for (int i = 2 * ntip + 2; i <= 2 * ntip + 2 * bc; i++) {

    ownFILL(owner, neighbour, clabel, -1, &bcode);
    GfaceFILL(faces, i, i + loop2, i + loop2 - 2, i - 2, &flabel, bcode,
              &oker[clabel], 4, &clabel);
    i++;
  };

  k = 1;
  while (k < j + nodes - 1) {

    for (int i = t2n + k * loop2; i < loop1 + k * loop2 - 2; i++) {
      ownFILL(owner, neighbour, clabel, -1, &bcode);
      GfaceFILL(faces, i, i + loop2, i + loop2 + 2, i + 2, &flabel, bcode,
                &oker[clabel], 4, &clabel);
      i++;
    };
    i1 = loop1 + k * loop2 - 2;
    comp = loop1 + (k - 1) * loop2;
    ownFILL(owner, neighbour, clabel, -1, &bcode);
    GfaceFILL(faces, i1, i1 + loop2, comp + loop2, comp, &flabel, bcode,
              &oker[clabel], 4, &clabel);

    for (int i = loop1 + (k - 1) * loop2;
         i < loop1 + (k - 1) * loop2 + 2 * ntip; i++) {

      ownFILL(owner, neighbour, clabel, -1, &bcode);
      GfaceFILL(faces, i, i + loop2, i + loop2 + 2, i + 2, &flabel, bcode,
                &oker[clabel], 4, &clabel);
      i++;
    };

    i1 = 2 * ntip + loop2 * k;
    comp = 2 * ntip + 2 * bc + 2 + loop2 * k;
    ownFILL(owner, neighbour, clabel, -1, &bcode);
    GfaceFILL(faces, i1, i1 + loop2, comp + loop2, comp, &flabel, bcode,
              &oker[clabel], 4, &clabel);

    for (int i = 2 * ntip + 2 * bc + 2 + loop2 * k; i < t2n + loop2 * k; i++) {

      ownFILL(owner, neighbour, clabel, -1, &bcode);
      GfaceFILL(faces, i, i + loop2, i + loop2 + 2, i + 2, &flabel, bcode,
                &oker[clabel], 4, &clabel);
      i++;
    };
    for (int i = 2 * ntip + 2 + loop1 + (k - 1) * loop2;
         i <= 2 * ntip + loop1 + (k - 1) * loop2 + 2 * bc; i++) {

      ownFILL(owner, neighbour, clabel, -1, &bcode);
      GfaceFILL(faces, i, i - 2, i - 2 + loop2, i + loop2, &flabel, bcode,
                &oker[clabel], 4, &clabel);
      i++;
    };

    for (int i = 2 * ntip + 4 + loop1 + (k - 1) * loop2 + 2 * bc;
         i <= 2 * ntip + 2 + loop1 + (k - 1) * loop2 + 4 * bc; i++) {

      ownFILL(owner, neighbour, clabel, -1, &bcode);
      GfaceFILL(faces, i, i + loop2, i + loop2 - 2, i - 2, &flabel, bcode,
                &oker[clabel], 4, &clabel);
      i++;
    };

    k++;
  };

  clabel = 0;
  BoundaryPatch p064;
  // Printing order is Front, Back, Aerofoil, Inlet, Outlet, Top1, Bottom1,
  // Top2, Bottom2

  // typedef enum { FrontAndBack, Aerofoil, Inlet, Outlet, Top1, Bottom1, Top2,
  // Bottom2} BoundaryPatch;

  /* startFace[Aerofoil] =  n_InternalFaces + nFaces[FrontAndBack];
   startFace[Inlet] =  n_InternalFaces + nFaces[FrontAndBack] +
   nFaces[Aerofoil]; startFace[Outlet] =   startFace[Inlet] + nFaces[Inlet];
   startFace[Top1] =  startFace[Outlet] + nFaces[Outlet];
   startFace[Bottom1] =  startFace[Top1] + nFaces[Top1];
   startFace[Top2] =  startFace[Bottom1] + nFaces[Bottom1];
   startFace[Bottom2] = startFace[Top2] + nFaces[Top2]; */

  int JMax = 4;

  char name11[8][21];
  char type11[8][21];
  char physicalType11[8][21];

  strcpy(name11[Aerofoil], "Aerofoil"), strcpy(name11[Inlet], "Inlet"),
      strcpy(name11[Outlet], "Outlet"),
      strcpy(name11[FrontAndBack], "FrontAndBack"),
      strcpy(name11[Top1], "Top1"), strcpy(name11[Bottom1], "Bottom1"),
      strcpy(name11[Top2], "Top2"), strcpy(name11[Bottom2], "Bottom2");

  strcpy(type11[Aerofoil], "wall"), strcpy(type11[Inlet], "patch"),
      strcpy(type11[Outlet], "patch"), strcpy(type11[FrontAndBack], "empty"),
      strcpy(type11[Top1], "patch"), strcpy(type11[Bottom1], "patch"),
      strcpy(type11[Top2], "patch"), strcpy(type11[Bottom2], "patch");

  strcpy(physicalType11[Aerofoil], "wall"),
      strcpy(physicalType11[Inlet], "inlet"),
      strcpy(physicalType11[Outlet], "outlet"),
      strcpy(physicalType11[FrontAndBack], "empty"),
      strcpy(physicalType11[Top1], "patch"),
      strcpy(physicalType11[Bottom1], "patch"),
      strcpy(physicalType11[Top2], "patch"),
      strcpy(physicalType11[Bottom2], "patch");

  // Printing Aerofoil

  printf("About to print to aerofoil\n");
  bcode = AEROFOIL;

  clabel = 0;
  for (int i = t2n; i < loop1 - 2; i += 2) {
    ownFILL(owner, neighbour, clabel, -1, &bcode);
    onlyFILL(faces, i, i + 2, i + 3, i + 1, &flabel, bcode);
    AerofoilF[flabel - startFace[Aerofoil]] = flabel;
    flabel++;
    clabel++;
  };

  i1 = loop1 - 2;
  comp = 0;
  ownFILL(owner, neighbour, clabel, -1, &bcode);
  onlyFILL(faces, i1, 0, 1, i1 + 1, &flabel, bcode);
  AerofoilF[flabel - startFace[Aerofoil]] = flabel;
  flabel++;
  clabel++;

  for (int i = 0; i < 2 * ntip; i += 2) {
    ownFILL(owner, neighbour, clabel, -1, &bcode);
    onlyFILL(faces, i, i + 2, i + 3, i + 1, &flabel, bcode);
    AerofoilF[flabel - startFace[Aerofoil]] = flabel;
    flabel++;
    clabel++;
  };
  i1 = 2 * ntip;
  comp = 2 * ntip + 2 * bc + 2;
  ownFILL(owner, neighbour, clabel, -1, &bcode);
  onlyFILL(faces, i1, comp, comp + 1, i1 + 1, &flabel, bcode);
  AerofoilF[flabel - startFace[Aerofoil]] = flabel;
  flabel++;
  clabel++;
  for (int i = 2 * ntip + 2 * bc + 2; i < t2n; i += 2) {

    ownFILL(owner, neighbour, clabel, -1, &bcode);
    onlyFILL(faces, i, i + 2, i + 3, i + 1, &flabel, bcode);
    AerofoilF[flabel - startFace[Aerofoil]] = flabel;
    flabel++;
    clabel++;
  };

  printf("\n Aerofoil Printed %d faces. nFaces[Aerofoil] = %d\n",
         flabel - startFace[Aerofoil], nFaces[Aerofoil]);
  // Printing Inlet
  bcode = INLET;

  k = j + nodes - 1;
  clabel = (j + nodes - 2) * (2 * ntip + 2 * bc);

  for (int i = t2n + k * loop2; i < loop1 + k * loop2 - 2; i += 2) {
    ownFILL(owner, neighbour, clabel, -1, &bcode);
    onlyFILL(faces, i, i + 1, i + 3, i + 2, &flabel, bcode);
    InletF[flabel - startFace[Inlet]] = flabel;
    flabel++;
    clabel++;
  }

  i1 = loop1 + k * loop2 - 2;
  comp = loop1 + (k - 1) * loop2;
  ownFILL(owner, neighbour, clabel, -1, &bcode);
  onlyFILL(faces, i1, i1 + 1, comp + 1, comp, &flabel, bcode);
  InletF[flabel - startFace[Inlet]] = flabel;
  flabel++;
  clabel++;

  for (int i = loop1 + (k - 1) * loop2;
       i < loop1 + (k - 1) * loop2 + 2 * ncrit + 2 * mp; i += 2) {
    ownFILL(owner, neighbour, clabel, -1, &bcode);
    onlyFILL(faces, i, i + 1, i + 3, i + 2, &flabel, bcode);
    InletF[flabel - startFace[Inlet]] = flabel;
    flabel++;
    clabel++;
  };

  printf("\n Inlet Printed %d faces. nFaces[Inlet] = %d\n",
         flabel - startFace[Inlet], nFaces[Inlet]);

  // Printing Outlet

  bcode = OUTLET;

  clabel = (j + nodes - 2) * (2 * ntip + 2 * bc) + 2 * ntip + 2 * bc - 1;

  k = j + nodes - 1;

  for (int i = 2 * ntip + 2 * bc + k * loop2; i > 2 * ntip + 2 * bc;
       i -= loop2) {

    ownFILL(owner, neighbour, clabel, -1, &bcode);
    onlyFILL(faces, i, i - loop2, i - loop2 + 1, i + 1, &flabel, bcode);
    OutletF[flabel - startFace[Outlet]] = flabel;
    flabel++;
    clabel = clabel - (2 * ntip + 2 * bc);
  };

  i1 = 2 * ntip + 2 * bc;
  comp = i1 + loop1;
  clabel = 2 * ntip + bc - 1;
  ownFILL(owner, neighbour, clabel, -1, &bcode);
  onlyFILL(faces, i1, comp, comp + 1, i1 + 1, &flabel, bcode);
  OutletF[flabel - startFace[Outlet]] = flabel;
  flabel++;

  k = j + nodes - 2;
  clabel += 2 * ntip + 2 * bc;
  for (int i = comp; i < comp + loop2 * k; i += loop2) {

    ownFILL(owner, neighbour, clabel, -1, &bcode);
    onlyFILL(faces, i, i + loop2, i + loop2 + 1, i + 1, &flabel, bcode);
    OutletF[flabel - startFace[Outlet]] = flabel;
    flabel++;
    clabel = clabel + (2 * ntip + 2 * bc);
  };

  printf("\n Outlet Printed %d faces. nFaces[Outlet] = %d\n",
         flabel - startFace[Outlet], nFaces[Outlet]);

  // Printing top1
  bcode = TOP;

  clabel = (j + nodes - 2) * (2 * ntip + 2 * bc) + mp + (mp * r_num) / r_dem +
           2 * ncrit;

  k = j + nodes - 2;

  for (int i = loop1 + k * loop2 + 2 * (mp + ncrit);
       i < loop1 + k * loop2 + 2 * ntip; i += 2) {

    ownFILL(owner, neighbour, clabel, -1, &bcode);
    onlyFILL(faces, i, i + 1, i + 3, i + 2, &flabel, bcode);
    Top1F[flabel - startFace[Top1]] = flabel;
    flabel++;
    clabel++;
  };

  printf("\n Top1 Printed %d faces. nFaces[Top1] = %d\n",
         flabel - startFace[Top1], nFaces[Top1]);

  // Printing bottom1
  bcode = BOTTOM;

  i1 = loop1 + k * loop2 + 2 * ntip + 2 * bc + 2;
  comp = loop1 + k * loop2 + 2 * ntip + 4 * bc + 4;

  ownFILL(owner, neighbour, clabel, -1, &bcode);
  onlyFILL(faces, i1, i1 + 1, comp + 1, comp, &flabel, bcode);
  Bottom1F[flabel - startFace[Bottom1]] = flabel;
  flabel++;
  clabel++;

  for (int i = comp; i < t2n + (k + 1) * loop2; i += 2) {
    ownFILL(owner, neighbour, clabel, -1, &bcode);
    onlyFILL(faces, i, i + 1, i + 3, i + 2, &flabel, bcode);
    Bottom1F[flabel - startFace[Bottom1]] = flabel;
    flabel++;
    clabel++;
  }

  printf("\n Bottom1 Printed %d faces. nFaces[Bottom1] = %d\n",
         flabel - startFace[Bottom1], nFaces[Bottom1]);

  // Printing Top2
  bcode = TOP2;

  for (int i = loop1 + k * loop2 + 2 * ntip;
       i < loop1 + k * loop2 + 2 * ntip + 2 * bc; i += 2) {

    ownFILL(owner, neighbour, clabel, -1, &bcode);
    onlyFILL(faces, i, i + 1, i + 3, i + 2, &flabel, bcode);
    Top2F[flabel - startFace[Top2]] = flabel;
    flabel++;
    clabel++;
  };

  printf("\n Top2 Printed %d faces. nFaces[Top2] = %d\n",
         flabel - startFace[Top2], nFaces[Top2]);

  // Printing Bottom2
  bcode = BOTTOM2;

  for (int i = i1; i < i1 + 2 * bc; i += 2) {
    ownFILL(owner, neighbour, clabel, -1, &bcode);
    onlyFILL(faces, i, i + 2, i + 3, i + 1, &flabel, bcode);
    Bottom2F[flabel - startFace[Bottom2]] = flabel;
    flabel++;
    clabel++;
  };

  printf("\n Bottom2 Printed %d faces. nFaces[Bottom2] = %d\n",
         flabel - startFace[Bottom2], nFaces[Bottom2]);

  printf("\nLast cell printed to is (4) %d\n", clabel - 1);

  cells[0] = fopen("cells", "a");
  cells[1] = fopen("cellcheck", "w");
  cells[2] = fopen("citago", "w");
  fprintf(cells[0], "\n\n");
  fprintf(cells[0], "%d\n(\n", sik);
  printf("Files opened successfully\n");

  boundar[0] = fopen("boundary", "a");
  boundar[1] = fopen("Truebound", "w");

  fprintf(boundar[0], "\n%d\n(\n", JMax);

  printf("\nAbout to print boundary\n");
  for (int J = 0; J < JMax; J++) {
    fprintf(boundar[0],
            "%s\n{\n\t type %s;\n \t physicalType %s;\n \t nFaces %d;\n \t "
            "startFace %d;\n } \n\n",
            name11[J], type11[J], physicalType11[J], nFaces[J], startFace[J]);
    fprintf(boundar[1], "%s\t%d\t%d\n", GetBoundaryName(J), nFaces[J],
            startFace[J]);
  };
  fprintf(boundar[0], ")\n\n// "
                      "********************************************************"
                      "***************** //");
  printf("Boundary printed correctly\n");

  for (int I5 = 0; I5 < sik; I5++) {
    fprintf(cells[0], "6(%d %d %d %d %d %d)\n", oker[I5].faces[0],
            oker[I5].faces[1], oker[I5].faces[2], oker[I5].faces[3],
            oker[I5].faces[4], oker[I5].faces[5]);

    fwrite(oker[I5].faces, sizeof(int), 6, cells[1]);
    fwrite(&oker[I5].label, sizeof(int), 1, cells[1]);
    // fprintf(cells[1], " %d %d %d %d %d %d\t%d\n", oker[I].faces[0],
    // oker[I].faces[1], oker[I].faces[2], oker[I].faces[3], oker[I].faces[4],
    // oker[I].faces[5], oker[I].label);

    fprintf(cells[2], "( %d %d %d %d %d %d )\t%d\n", oker[I5].faces[0],
            oker[I5].faces[1], oker[I5].faces[2], oker[I5].faces[3],
            oker[I5].faces[4], oker[I5].faces[5], oker[I5].label);
  };
  printf("Cells printed correctly\n");
  fprintf(cells[0], "\n)\n");
  fprintf(faces[0], "\n)\n");
  fprintf(neighbour[0], "\n)\n");
  fprintf(owner[0], "\n)\n");
  // End Here

  printf("\n t val of corner = %.15Lf", foil->associatedT[ncrit + mp]);
  printf("\n t val of bottom left corner = %.15Lf",
         foil->associatedT[ncrit + (int)(2 * mp / 3)]);
  printf("\n last element of foil count is %.15Lf",
         foil->associatedT[foil->count - ncrit - (int)(2 * mp / 3)]);
  printf("\n ncrit = %d", ncrit);
  printf("\n n tip = %d", array.size - 1);
  printf("\n cparfinal = %.15Lf", cpar);
  printf("\n tracecount = %d", tracecount);
  printf("\n j = %d", j);
  printf("\n nodes = %d", nodes);
  printf("\n loop1 = %d", loop1);
  printf("\n loop2 = %d", loop2);
  printf("\n CELLmod = %d", 2 * ntip + 2 * bc);
  printf("\n Total points = loop1 + (j + nodes - 1) loop2 = %d",
         loop1 + (j + nodes - 1) * loop2);
  printf("\n t2n =  %d",
         2 * ntip + 2 * backtrace + 2 * ntip - 2 * ncrit -
             2 * ((int)((((long double)(r_num)) / ((long double)(r_dem))) *
                        (mp))));
  long double tmpt = foil->associatedT[ncrit + (int)(2 * mp / 3)];
  printf("\n t(x,y) val of bottom left corner t2 = %.15Lf\t%.15Lf\n",
         tmpt + f3(tmpt, a) * sinl(theta1(tmpt, a)),
         f1(tmpt, a) - f3(tmpt, a) * cosl(theta1(tmpt, a)));

  printf("\n Total faces calculated: %d, Total faces: %d",
         2 * (2 * ntip + 2 * bc) * (j + nodes - 1) + (2 * ntip + bc) +
             (2 * ntip + 2 * bc) * (j + nodes - 1) +
             (2 * ntip + 1 + 2 * bc) * (j + nodes - 1),
         flabel);

  printf("\n n value of t2, saved as t2n %d", t2n);
  printf("\n n value of t1 %d", 2 * ncrit + 2 * mp);
  printf("\n t(x,y) val of corner t1 = %.15Lf\t%.15Lf\n",
         foil->associatedT[ncrit + mp] -
             f3(foil->associatedT[ncrit + mp], a) *
                 sinl(theta1(foil->associatedT[ncrit + mp], a)),
         f1(foil->associatedT[ncrit + mp], a) +
             f3(foil->associatedT[ncrit + mp], a) *
                 cosl(theta1(foil->associatedT[ncrit + mp], a)));
  printf("\nTotal cell count is allegedly %d\n", sik);
  printf("\np1 = %.7LF, p2 = %.7Lf\n",
         (foil->associatedT[ntip - 1] - a.c * a.p) /
             (foil->associatedT[ntip] - a.c * a.p),
         1 - (foil->associatedT[ntip - 1] - a.c * a.p) /
                 (foil->associatedT[ntip] - a.c * a.p));

  fclose(faces[0]);
  fclose(faces[1]);
  fclose(faces[2]);
  fclose(cells[0]);
  fclose(cells[1]);
  fclose(cells[2]);
  free(oker);
};

void foilcompute(FILE *pointoutput[3], foilpoint *foil, aerofoil a,
                 Tarray *array, double tmax, int backtrace, trace *trace1,
                 long double xR) {
  pointoutput[0] = fopen("figaro1.txt", "w");
  pointoutput[1] = fopen("figarosafe.txt", "w");
  pointoutput[2] = fopen("points", "a");
  pointoutput[3] = fopen("pointcheck", "w");

  if ((pointoutput[0] == NULL || pointoutput[1] == NULL)) {
    printf("\nERROR OPENING FILE\n");
  };
  fprintf(pointoutput[2], "\n\n(\n");
  foil->count = 2 * array->size - 2;
  printf("\nfoilcount = %d\n", foil->count);
  trace1->size = backtrace;
  trace1->x = (long double *)malloc((trace1->size) * sizeof(long double));

  trace1->label = (int *)malloc((trace1->size) * sizeof(int));

  foil->x = (long double *)malloc((foil->count) * sizeof(long double));

  foil->y = (long double *)malloc((foil->count) * sizeof(long double));

  foil->slope = (long double *)malloc((foil->count) * sizeof(long double));

  foil->associatedT =
      (long double *)malloc((foil->count) * sizeof(long double));
  foil->label = (int *)malloc((foil->count) * sizeof(int));

  foil->section = (int *)malloc((foil->count) * sizeof(int));
  int n = 0;
  int cable = 0;
  long double z0 = 0.0L;

  while (*(array->T + n) < a.p * a.c) {
    long double tmpt;
    tmpt = *(array->T + n);
    foil->x[n] = tmpt - f3(tmpt, a) * sin(theta1(tmpt, a));
    foil->y[n] = f1(tmpt, a) + f3(tmpt, a) * cos(theta1(tmpt, a));
    foil->label[n] = n;
    foil->section[n] = 1;
    foil->associatedT[n] = tmpt;
    if (n == 0) {
      foil->slope[0] = (-1.0 / dydx1(tmpt + 0.000000000001, a));
    } else {
      foil->slope[n] = (-1.0 / dydx1(tmpt, a));
    }
    fprintf(pointoutput[0], "%.15Lf\t%.15Lf\n", foil->x[n], foil->y[n]);
    //   printf("t(%d) = %.15Lf\t%.15Lf\t%.15Lf\n", foil->label[n],
    //   foil->associatedT[n], foil->x[n], foil->y[n]);
    fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[n],
            foil->y[n], 0.0L, 2 * foil->label[n]);

    fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[n],
            foil->y[n], delz, 2 * foil->label[n] + 1);

    WriteToPoint(&foil->x[n], &foil->y[n], &(long double){0.0L},
                 &(int){2 * foil->label[n]}, pointoutput[3]);
    WriteToPoint(&foil->x[n], &foil->y[n], &(long double){delz},
                 &(int){2 * foil->label[n] + 1}, pointoutput[3]);

    fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[n], foil->y[n],
            0.0L);
    fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[n], foil->y[n],
            delz);

    n++;
  };

  while (n < array->size) {
    long double tmpt;
    tmpt = *(array->T + n);
    foil->x[n] = tmpt - f3(tmpt, a) * sin(theta2(tmpt, a));
    foil->y[n] = f2(tmpt, a) + f3(tmpt, a) * cos(theta2(tmpt, a));
    foil->section[n] = 3;
    foil->label[n] = n;
    foil->associatedT[n] = tmpt;

    foil->slope[n] = (-1.0 / dydx3(tmpt, a));

    fprintf(pointoutput[0], "%.15Lf\t%.15Lf\n", foil->x[n], foil->y[n]);
    //   printf("t(%d) = %.15Lf\t%.15Lf\t%.15Lf\n", foil->label[n],
    //   foil->associatedT[n], foil->x[n], foil->y[n]);

    fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[n],
            foil->y[n], 0.0L, 2 * foil->label[n]);
    fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[n],
            foil->y[n], delz, 2 * foil->label[n] + 1);

    WriteToPoint(&foil->x[n], &foil->y[n], &(long double){0.0L},
                 &(int){2 * foil->label[n]}, pointoutput[3]);
    WriteToPoint(&foil->x[n], &foil->y[n], &(long double){delz},
                 &(int){2 * foil->label[n] + 1}, pointoutput[3]);

    fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[n], foil->y[n],
            0.0L);
    fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[n], foil->y[n],
            delz);

    if (n == array->size - 1) {
      long double deltax = (xR - *(array->T + n)) / ((long double)(backtrace));
      for (int i = 1; i <= backtrace; i++) {
        trace1->x[i - 1] = ((long double)(i)) * deltax + foil->x[n];
        trace1->y = foil->y[n];
        trace1->label[i - 1] = n + i;
        fprintf(pointoutput[0], "%.15Lf\t%.15Lf\n", trace1->x[i - 1],
                trace1->y);

        fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", trace1->x[i - 1],
                trace1->y, 0.0L, 2 * trace1->label[i - 1]);
        fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", trace1->x[i - 1],
                trace1->y, delz, 2 * trace1->label[i - 1] + 1);

        WriteToPoint(&trace1->x[i - 1], &trace1->y, &(long double){0.0L},
                     &(int){2 * trace1->label[i - 1]}, pointoutput[3]);
        WriteToPoint(&trace1->x[i - 1], &trace1->y, &(long double){delz},
                     &(int){2 * trace1->label[i - 1] + 1}, pointoutput[3]);

        fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", trace1->x[i - 1],
                trace1->y, 0.0L);
        fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", trace1->x[i - 1],
                trace1->y, delz);
      };

      printf("\n tip %.15Lf occurs at %d\n", *(array->T + n), n);
    };
    n++;
  };
  int k = n;
  n = array->size - 2;

  while (*(array->T + n) > a.p * a.c) {

    long double tmpt;

    tmpt = *(array->T + n);

    foil->x[k] = tmpt + f3(tmpt, a) * sin(theta2(tmpt, a));
    foil->y[k] = f2(tmpt, a) - f3(tmpt, a) * cos(theta2(tmpt, a));
    foil->section[k] = 4;
    foil->label[k] = k + backtrace;

    foil->slope[k] = (-1.0 / dydx4(tmpt, a));
    foil->associatedT[k] = tmpt;
    fprintf(pointoutput[0], "%.15Lf\t%.15Lf\n", foil->x[k], foil->y[k]);
    //   printf("t(%d) = %.15Lf\t%.15Lf\t%.15Lf\n", foil->label[k],
    //   foil->associatedT[k], foil->x[k], foil->y[k]);

    fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[k],
            foil->y[k], 0.0L, 2 * foil->label[k]);
    fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[k],
            foil->y[k], delz, 2 * foil->label[k] + 1);

    WriteToPoint(&foil->x[k], &foil->y[k], &(long double){0.0L},
                 &(int){2 * foil->label[k]}, pointoutput[3]);
    WriteToPoint(&foil->x[k], &foil->y[k], &(long double){delz},
                 &(int){2 * foil->label[k] + 1}, pointoutput[3]);

    fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[k], foil->y[k],
            0.0L);
    fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[k], foil->y[k],
            delz);

    n--;
    k++;
  };

  while (n > 0) {
    long double tmpt;
    tmpt = *(array->T + n);
    foil->x[k] = tmpt + f3(tmpt, a) * sin(theta1(tmpt, a));
    foil->y[k] = f1(tmpt, a) - f3(tmpt, a) * cos(theta1(tmpt, a));
    foil->section[k] = 2;
    foil->label[k] = k + backtrace;

    foil->associatedT[k] = tmpt;

    foil->slope[k] = (-1.0 / dydx2(tmpt, a));

    fprintf(pointoutput[0], "%.15Lf\t%.15Lf\n", foil->x[k], foil->y[k]);
    //   printf("t(%d) = %.15Lf\t%.15Lf\t%.15Lf\n", foil->label[k],
    //   foil->associatedT[k], foil->x[k], foil->y[k]);

    fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[k],
            foil->y[k], 0.0L, 2 * foil->label[k]);
    fprintf(pointoutput[1], "%.15Lf\t%.15Lf\t%.4Lf\t%d\n", foil->x[k],
            foil->y[k], delz, 2 * foil->label[k] + 1);

    WriteToPoint(&foil->x[k], &foil->y[k], &(long double){0.0L},
                 &(int){2 * foil->label[k]}, pointoutput[3]);
    WriteToPoint(&foil->x[k], &foil->y[k], &(long double){delz},
                 &(int){2 * foil->label[k] + 1}, pointoutput[3]);

    fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[k], foil->y[k],
            0.0L);
    fprintf(pointoutput[2], "(%.15Lf %.15Lf %.4Lf )\n", foil->x[k], foil->y[k],
            delz);
    n--;
    k++;
  };
  printf("Array size was %d\n", array->size);
  printf("Allocated point size was %d", foil->count);
  fclose(pointoutput[0]);
  fclose(pointoutput[1]);
  fclose(pointoutput[2]);
  fclose(pointoutput[3]);
};

void Tvalues(int mp, long double r, long double t[4], Tarray *array,
             int *defaulter, int gp) {
  long double A, B, C;
  long double A1, B1;

  C = t[0];
  printf("C = %.15Lf\n", C);
  if (t[2] > t[0]) {
    B = log((t[2] - t[0]) / (t[1] - t[0])) * (1 / log(r));
    defaulter[0] = 0;
  } else {
    B = log((t[0] - t[2]) / (t[1] - t[0])) * (1 / log(r));
    defaulter[0] = 1;
  };
  printf("B = %.15Lf\n", B);
  A = (t[1] - t[0]) / powl(((long double)(mp)), B);

  A1 = A * powl((mp + 1), B) + t[0] - t[1];

  B1 = (logl(fabsl((t[3] - t[1]) / A1)) / (logl(gp)));

  long double temp = 0;
  int n = 0;
  printf("A = %.15Lf\n", A);
  printf("\nThe size is %d\n", ((int)(pow((t[3] - C) / A, (1.0 / B)))) + 2);

  array->size = ((int)(pow((t[3] - C) / A, (1.0 / B)))) + 2;

  array->size = mp + gp + 1;

  printf("%d\n", array->size);

  array->T = (long double *)malloc((array->size) * sizeof(long double));

  while (n <= mp) {
    temp = A * pow((long double)n, B) + C;
    array->T[n] = temp;
    //   printf("t(%d) = %.15Lf\n", n, array->T[n]);
    n++;
  };
  while (n < array->size && n > mp) {
    int m11 = n - mp;
    temp = A1 * pow((long double)m11, B1) + t[1];
    array->T[n] = temp;
    //    printf("t(%d) = %.15Lf\n", n, array->T[n]);
    n++;
  };

  printf(" t  %d ] = %.15Lf\n", n - 1, t[3]);
  printf("\nThe size %d\n", array->size);
  printf("\n ntip at %d is %.15Lf", array->size - 1, array->T[array->size - 1]);
  array->T[array->size - 1] = t[3];

  //  printf("tf(%d) = %.15Lf\n", array->size - 1, array->T[array->size - 1]);
  //  printf("size still is\t\n%d", array->size);
};

long double arcelength(long double t1, long double t2, int whichfunction) {
  // x1 = t - f3(t) * sin (theta1(t) );
  // y1 = f1(t) + f3(t) cos (theta1(t));
  // x1' = 1 - f3prime(t) sin (theta1(t)) - f3(t) theta1prime(t) cos(theta1(t))
  // y1' = f1prime(t) + f3prime(t) cos (theta1(t)) - f3(t) theta1prime(t) sin
  // (theta1(t))

  // Integrand = sqrt(x1'^2 + y1'^2)
  //
  // w_i = 2 / ( (1 - x_i)^2 [P'n(x_i)]^2 )
  // P_n(x) = (1/(n! 2^n)) d^n/dx^n (x^2 - 1)^n
  FILE *generatedc;
  generatedc = fopen("generated.c", "w");
  int n = 2;

  fprintf(generatedc, "#include <stdio.h>\n");
  fprintf(generatedc, "#include <math.h>\n");

  switch (whichfunction) {
  case 1:

    break;
  default:
    break;
  };

  return 1.0L;
};

void arcadjust(Tarray *array, aerofoil a, int *ncrit, int *defaulter) {
  long double t0 = array->T[0];
  long double y0 =
      f1(array->T[0], a) + f3(array->T[0], a) * cosl(theta1(array->T[0], a));
  long double y1 =
      f1(array->T[1], a) + f3(array->T[1], a) * cosl(theta1(array->T[1], a));
  long double x0 =
      array->T[0] - f3(array->T[0], a) * sin(theta1(array->T[0], a));
  long double x1 =
      array->T[1] - f3(array->T[1], a) * sin(theta1(array->T[1], a));
  long double d1, d2;
  d1 = sqrtl(powl((y1 - y0), 2) + powl((x1 - x0), 2));
  d2 = sqrtl(powl(y0, 2) + powl(x0, 2));
  printf("/////// arcadjust data ///////\n");
  printf("(x0, y0) = (%.7Lf, %.7Lf)\n", x0, y0);
  printf("(x1, y1) = (%.7Lf, %.7Lf)\n", x1, y1);
  printf("(d1, d2) = (%.7Lf, %.7Lf)\n", d1, d2);

  int fefOG = ((int)(d2 / d1)) + 1;
  printf("fefOG = %d\n", fefOG);

  if (fefOG > 1 && (*defaulter == 0)) {
    long double A = array->T[0] / (powl((long double)(fefOG), 3));

    array->T = realloc(array->T, (array->size + fefOG) * sizeof(long double));

    memmove(array->T + fefOG, array->T, (array->size) * sizeof(long double));
    //////////////////////////////////////////////////////////
    for (int i = 0; i < fefOG; i++) {
      array->T[i] = (A * powl(((long double)(i)), 3));
    };
    array->size = array->size + fefOG;
  };
  if (fefOG == 1 || (defaulter == 0)) {
    array->T = realloc(array->T, (array->size + fefOG) * sizeof(long double));
    memmove(array->T + fefOG, array->T, (array->size) * sizeof(long double));
    array->T[0] = 0;
    array->size += 1;
  };

  *ncrit = fefOG;
  printf("\nPay attention: fefOG = %d \t whilst ncrit =%d ", fefOG, *ncrit);
  printf("\nOld Tarray size: %d\n", array->size - fefOG);
  printf("\nNew size = %d", array->size);
}
#define PREC 256

void PrintData(MeshFineness *Mesh1, aerofoil *aero) {
  printf("xL = %.15Lf \n", Mesh1->xL);
  printf("xR = %.15Lf \n", Mesh1->xR);
  printf("yB = %.15Lf \n", Mesh1->yB);
  printf("yU = %.15Lf \n", Mesh1->yU);
  printf("m =  %.15Lf \n", aero->m);
  printf("P =  %.15Lf \n", aero->p);
  printf("T =  %.15Lf \n", aero->T);
  printf("c =  %.15Lf \n", aero->c);
  printf("mp = %d  \n", Mesh1->mp);
  printf("r_num = %d  \n", Mesh1->r_num);
  printf("r_dem = %d  \n", Mesh1->r_dem);
  printf("gp =    %d    \n", Mesh1->gp);
  printf("cpar =  %.6Lf  \n", Mesh1->Cpar);
  printf("tpar =  %d \n", Mesh1->Tpar);
};

int readMesh2(FILE *inputfile, MeshFineness *Mesh1, aerofoil *aero) {

  int j = 0;
  j = fscanf(inputfile, "%Lf\n", &Mesh1->xL);
  j += fscanf(inputfile, "%Lf\n", &Mesh1->xR);
  j += fscanf(inputfile, "%Lf\n", &Mesh1->yB);
  j += fscanf(inputfile, "%Lf\n", &Mesh1->yU);
  j += fscanf(inputfile, "%Lf\n", &aero->m);
  j += fscanf(inputfile, "%Lf\n", &aero->p);
  j += fscanf(inputfile, "%Lf\n", &aero->T);
  j += fscanf(inputfile, "%Lf\n", &aero->c);
  j += fscanf(inputfile, "%d\n", &Mesh1->mp);
  j += fscanf(inputfile, "%d\n", &Mesh1->r_num);
  j += fscanf(inputfile, "%d\n", &Mesh1->r_dem);
  j += fscanf(inputfile, "%d\n", &Mesh1->gp);
  j += fscanf(inputfile, "%Lf\n", &Mesh1->Cpar);
  j += fscanf(inputfile, "%d\n", &Mesh1->Tpar);
  PrintData(Mesh1, aero);
  return j;
}

int ReadMesh1(FILE *inputfile, MeshFineness *Mesh1, aerofoil *aero) {
  int j = 0;
  j += fscanf(inputfile, "xL = %Lf\n", &Mesh1->xL);
  j += fscanf(inputfile, "xR = %Lf\n", &Mesh1->xR);
  j += fscanf(inputfile, "yB = %Lf\n", &Mesh1->yB);
  j += fscanf(inputfile, "yU = %Lf\n", &Mesh1->yU);
  j += fscanf(inputfile, "m = %Lf\n", &aero->m);
  j += fscanf(inputfile, "P = %Lf\n", &aero->p);
  j += fscanf(inputfile, "T = %Lf\n", &aero->T);
  j += fscanf(inputfile, "c = %Lf\n", &aero->c);
  j += fscanf(inputfile, "mp = %d\n", &Mesh1->mp);
  j += fscanf(inputfile, "r_num = %d\n", &Mesh1->r_num);
  j += fscanf(inputfile, "r_dem = %d\n", &Mesh1->r_dem);
  j += fscanf(inputfile, "gp = %d\n", &Mesh1->gp);
  j += fscanf(inputfile, "cpar = %Lf\n", &Mesh1->Cpar);
  j += fscanf(inputfile, "tpar = %d\n", &Mesh1->Tpar);
  printf("\nHere's what ReadMesh1 Found: \n");
  PrintData(Mesh1, aero);
  return j;
};

void infoprint() {
  printf("xL = ...\n");
  printf("xR = ...\n");
  printf("yB = ...\n");
  printf("yU = ...\n");
  printf("m = ...\n");
  printf("P = ...\n");
  printf("T = ...\n");
  printf("c = ...\n");
  printf("mp = ...\n");
  printf("r_num = ...\n");
  printf("r_dem = ...\n");
  printf("gp = ...\n");
  printf("cpar = ...\n");
  printf("tpar = ...\n");
}
void descriptionprint() {
  printf("\n\n Where xL is the leftmost coordinate of the mesh, xR is the "
         "rightmost coordiante of the mesh, yU is the uppermost coordinate of "
         "the mesh\n");
  printf("m is the camber, P is the position of the maximum camber, T is the "
         "thickness, c is the chord length\n");
  printf(
      "mp is the number of points between the horizontal median point and the "
      "point which connects the top left point of the mesh\n mp x r_num/r_dem "
      "is the number of points between the horizontal median point and the "
      "point which connects the bottom left point of the mesh\n, gp is the "
      "number of points between the mp-th point and the tip of the airfoil\n, "
      "cpar is the coefficient which dicates spacing of points\n, and tpar is "
      "the coefficient which dictates the density of the tail.\n");
};

int main(int argc, char *argv[]) {
  setlocale(LC_ALL, "");

  long double xL, xR, yU, yB;
  int mp;
  long double r = (4.0 / 5.0);
  int r_num = 4;
  int r_dem = 5;
  int precision;
  long double t0, t1, t2, t3;
  int kok;
  aerofoil aero;
  MeshFineness Mesh1;

  long double a_vals[] = {0.2969, -.1260, -.3516, .2843, -.1015};

  FILE *inputfile;
  unsigned int FileProvided = 0;
  unsigned int FileType = 0x9283A942;

  switch (argc) {
  case 1:
    FileProvided = 0;
    break;
  case 2:
  badcase0:
    printf("Either 0 or 2 arguments is allowed. Format is as follows\n./b "
           "[File type] [File name]\n");
    printf("Or, you can use ./b by itself to enter values manually\n");
    if (FileType != 0x9283A942) {
      printf("You've entered a File type of %d. Valid options are 0 or 1\n",
             FileType);
    };
    printf("File type 1 example:\n");
    infoprint();
    printf("File type 0 equivalent example:\n");
    printf("-5\n10\n-10\n10\n0.02\n0.4\n0.11\n1\n12\n2\n3\n60\n0.123\n5\n");

    return 4;
    break;
  case 3:
    inputfile = fopen(argv[2], "r");
    printf("Reading file ");
    printf("%s\n", argv[2]);
    FileProvided = 1;
    FileType = atoi(argv[1]);
    break;
  default:
    printf("Too many arguments. Either 0 or 2 arguments is allowed. Format is "
           "either:\n");
    printf("./b\n");
    printf("or\n");
    printf("./b [File type] [File name]\n");
    printf("File type 1 example:\n");
    infoprint();
    printf("File type 0 equivalent example:\n");
    printf("-5\n10\n-10\n10\n0.02\n0.4\n0.11\n1\n12\n2\n3\n60\n0.123\n5\n");
  badcase1:

    return 4;
  };

  memcpy(&aero.a, a_vals, sizeof(a_vals));

  printf("%s", "Please Enter the left most part of the mesh, xL:\n");
  Mesh1.xL = -5;
  printf("%s", "Please Enter the right most part of the mesh, xR:\n");
  Mesh1.xR = 10;
  printf("%s", "Please Enter the Upper part of the mesh, yU:\n");
  Mesh1.yU = 10;
  printf("%s", "Please Enter the Bottom most part of the mesh, yB:\n");
  Mesh1.yB = -10;
  xL = Mesh1.xL;
  xR = Mesh1.xR;
  yU = Mesh1.yU;
  yB = Mesh1.yB;

  int op59 = 0;

  if (FileProvided == 1 && FileType == 0) {
    op59 = readMesh2(inputfile, &Mesh1, &aero);
  };
  if (FileProvided == 1 && FileType == 1) {
    op59 = ReadMesh1(inputfile, &Mesh1, &aero);
  };
  printf("op59 = %d\n", op59);

  if (FileProvided == 1 && op59 != 14) {
    printf("File was provided, however op59 = %d\n", op59);
    // goto badcase1;
  };

  printf("%s", "Enter m, then p, then T, then c\n");

  if ((FileProvided == 1 && op59 != 14) || (FileProvided == 0)) {
    printf("%s", "Enter m, then p, then T, then c\n");
    printf("m = ");
    scanf("%Lf", &aero.m);
    printf("p = ");
    scanf("%Lf", &aero.p);
    printf("T = ");
    scanf("%Lf", &aero.T);
    printf("c = ");
    scanf("%Lf", &aero.c);
  }

  printf("m = %.2Lf\np = %.2Lf\nT = %.2Lf\n", aero.m, aero.p, aero.T);

  precision = 15;

  aero.kk = aero.m / (aero.p * aero.p);
  aero.jj = aero.m / ((1 - aero.p) * (1 - aero.p));

  int j = 0;
  long double prevtval = 0.000000000000000001;
  long double correctdigit = 0;
  long double testval = 0;
  long double minval = 1;
  long double fixdigit = 0;
  long double deroundee;
  int refiner = 0;

  while (j <= precision + 3) {

    for (int i = 0; i <= 9; i++) {
      long double tempval[10];
      testval = prevtval + ((long double)(i)) / pow(10.0, (long double)(j + 2));
      tempval[i] = Zeroroot(testval, aero);
      if (i == 0) {
        minval = tempval[i];
        correctdigit = ((long double)i) / pow(10.0, (long double)(j + 2));
        fixdigit = correctdigit;
      } else {
        if (fabsl(tempval[i]) < fabsl(minval)) {
          minval = tempval[i];
          correctdigit = ((long double)i) / pow(10.0, (long double)(j + 2));
          fixdigit = correctdigit;
        }
      }
    }
    if (fixdigit <= 1e-17) {
      goto zerocase;
    }
  justforfun:
    for (int i1 = -9; i1 <= 9; i1++) {
      testval = prevtval + fixdigit +
                ((long double)i1) / pow(10.0, (long double)(j + 3 + refiner));
      long double temp2val[19];
      temp2val[i1 + 9] = Zeroroot(testval, aero);
      if (i1 == -9) {
        deroundee = fixdigit + ((long double)i1) /
                                   pow(10.0, (long double)(j + 3 + refiner));
        minval = temp2val[i1 + 9];

      } else {
        if (fabsl(temp2val[i1 + 9]) < fabsl(minval)) {
          deroundee = fixdigit + ((long double)i1) /
                                     pow(10.0, (long double)(j + 3 + refiner));
          minval = temp2val[i1 + 9];
        }
        if ((i1 == 9) && (fabsl(deroundee - fixdigit) <= 1e-17L) &&
            (j + 3 + refiner < 16)) {
          refiner++;
          goto justforfun;
        }
      }
    }

    if (deroundee < fixdigit) {
      correctdigit = correctdigit - (1.0 / pow(10.0, (long double)(j + 2)));
    }
  zerocase:
    prevtval = prevtval + correctdigit;

    j++;
  }
  t0 = prevtval;
  printf("%s %.17Lf\n", "t0 = ", t0);

  j = 0;
t1fix:
  prevtval = 0.0000000000000000001;
  correctdigit = 0;
  testval = 0;
  minval = 0;
  fixdigit = 0;
  deroundee = 0;
  refiner = 0;
  while (j <= precision) {

    for (int i = 0; i <= 9; i++) {
      long double tempval[10];
      testval = prevtval + ((long double)i) / pow(10.0, (long double)(j + 2));
      tempval[i] = Firstroot(testval, aero, xL, yU);
      if (i == 0) {
        minval = tempval[i];
        correctdigit = ((long double)i) / pow(10.0, (long double)(j + 2));
        fixdigit = correctdigit;
      } else {
        if (fabsl(tempval[i]) < fabsl(minval)) {
          minval = tempval[i];
          correctdigit = ((long double)i) / pow(10.0, (long double)(j + 2));
          fixdigit = correctdigit;
        }
      }
    }
    if (fixdigit <= 1e-17L) {
      goto zerocase1;
    }
  justforfun1:
    for (int i1 = -9; i1 <= 9; i1++) {
      testval = prevtval + fixdigit +
                ((long double)i1) / pow(10.0, (long double)(j + 3 + refiner));

      long double temp2val[19];
      temp2val[i1 + 9] = Firstroot(testval, aero, xL, yU);
      if (i1 == -9) {
        minval = temp2val[i1 + 9];
        deroundee = fixdigit + ((long double)i1) /
                                   pow(10.0, (long double)(j + 3 + refiner));

      } else {
        if (fabsl(temp2val[i1 + 9]) < fabsl(minval)) {
          minval = temp2val[i1 + 9];
          deroundee = fixdigit + ((long double)i1) /
                                     pow(10.0, (long double)(j + 3 + refiner));
        }
      }
      if ((i1 == 9) && (fabsl(deroundee - fixdigit) <= 1e-17L) &&
          (refiner < 16)) {
        refiner++;
        goto justforfun1;
      }
    }

    if (deroundee < fixdigit) {
      correctdigit = correctdigit - (1.0 / pow(10.0, (long double)(j + 2)));
    }
  zerocase1:
    prevtval = prevtval + correctdigit;

    j++;
  }
  j = 0;
  t1 = prevtval;
  printf("%s %.18Lf\n", "t1 = ", t1);
  if (fabsl(Firstroot(t1, aero, xL, yU)) > 0.01) {
    int ii = 0;
    printf("First %.15Lf !=0. (2 to increment exponent, 1 to decrement "
           "exponent, 0 to ignore",
           Secondroot(t2, aero, xL, yB));
    scanf("%d", &ii);
    if (ii == 1) {
      j--;
      goto t1fix;
    };
    if (ii == 2) {
      j++;
      goto t1fix;
    };
  };

  j = 0;
t2fix:
  prevtval = 0.000000000000000001;
  refiner = 0;
  correctdigit = 0;
  testval = 0;
  minval = 0;
  fixdigit = 0;
  deroundee = 0;
  while (j <= precision + 3) {

    for (int i = 0; i <= 9; i++) {
      long double tempval[10];

      testval = prevtval + ((long double)i) / pow(10.0, (long double)(j + 2));
      tempval[i] = Secondroot(testval, aero, xL, yB);
      if (i == 0) {
        minval = tempval[i];
        correctdigit = ((long double)i) / pow(10.0, (long double)(j + 2));
        fixdigit = correctdigit;
      } else {
        if (fabsl(tempval[i]) < fabsl(minval)) {
          minval = tempval[i];
          correctdigit = ((long double)i) / pow(10.0, (long double)(j + 2));
          fixdigit = correctdigit;
        }
      }
    }
    if (fixdigit <= 1e-17) {
      goto zerocase2;
    }
  justforfun2:
    for (int i1 = -9; i1 <= 9; i1++) {
      testval = prevtval + fixdigit +
                ((long double)i1) / pow(10.0, (long double)(j + 3 + refiner));
      long double temp2val[19];
      temp2val[i1 + 9] = Secondroot(testval, aero, xL, yB);
      if (i1 == -9) {
        deroundee = fixdigit + ((long double)i1) /
                                   pow(10.0, (long double)(j + 3 + refiner));
        minval = temp2val[i1 + 9];

      } else {
        if (fabsl(temp2val[i1 + 9]) < fabsl(minval)) {
          deroundee = fixdigit + ((long double)i1) /
                                     pow(10.0, (long double)(j + 3 + refiner));
          minval = temp2val[i1 + 9];
        }
      }
      if ((i1 == 9) && (fabsl(deroundee - fixdigit) <= 1e-17L) &&
          (refiner < 16)) {
        refiner++;
        goto justforfun2;
      }
    }
    if (deroundee < fixdigit) {
      correctdigit = correctdigit - (1.0 / pow(10.0, (long double)(j + 2)));
    }
  zerocase2:
    prevtval = prevtval + correctdigit;

    j++;
  }
  t2 = prevtval;
  j = 0;
  prevtval = 1.0000000000000;
  correctdigit = 0;
  testval = 0;
  minval = 0;
  fixdigit = 0;
  deroundee = 0;
  refiner = 0;
  printf("%s %.15Lf\n", "t2 = ", t2);
  if (fabsl(Secondroot(t2, aero, xL, yB)) > 0.01) {
    int ii = 0;
    printf("Secondroot %.15Lf !=0. (2 to increment exponent, 1 to decrement "
           "exponent, 0 to ignore",
           Secondroot(t2, aero, xL, yB));
    scanf("%d", &ii);
    if (ii == 1) {
      j--;
      goto t2fix;
    };
    if (ii == 2) {
      j++;
      goto t2fix;
    };
  };

  printf("%s %.15Lf\n", "t2 = ", t2);

  while (j <= precision) {

    for (int i = 0; i <= 9; i++) {
      long double tempval[10];
      testval = prevtval + ((long double)i) / pow(10.0, (long double)(j));
      tempval[i] = Thirdroot(testval, aero);
      //      printf("Thirdroot ( %.15Lf ) = %.15Lf\n", testval, tempval[i]);
      if (i == 0) {
        minval = tempval[i];
        correctdigit = ((long double)i) / pow(10.0, (long double)(j));
        fixdigit = correctdigit;
      } else {
        if (fabsl(tempval[i]) < fabsl(minval)) {
          minval = tempval[i];
          correctdigit = ((long double)i) / pow(10.0, (long double)(j));
          fixdigit = correctdigit;
        }
      }
    }
    if (fixdigit <= 1e-17L) {
      goto zerocase3;
    }
  justforfun3:
    for (int i1 = -9; i1 <= 9; i1++) {
      testval = prevtval + fixdigit +
                ((long double)i1) / pow(10.0, (long double)(j + 1 + refiner));

      long double temp2val[19];
      temp2val[i1 + 9] = Thirdroot(testval, aero);
      if (i1 == -9) {
        deroundee = fixdigit + ((long double)i1) /
                                   pow(10.0, (long double)(j + 1 + refiner));
        minval = temp2val[i1 + 9];

      } else {
        if (fabsl(temp2val[i1 + 9]) < fabsl(minval)) {
          deroundee = fixdigit + ((long double)i1) /
                                     pow(10.0, (long double)(j + 1 + refiner));
          minval = temp2val[i1 + 9];
        }
      }
      if ((i1 == 9) && (fabsl(deroundee - fixdigit) <= 1e-17L) &&
          (refiner < 16)) {
        printf(
            "The digit was found to be 0. Heres what it is %.15Lf = %.15Lf\n",
            deroundee, fixdigit);
      }
      if ((i1 == 9) && (fabsl(deroundee - fixdigit) <= 1e-17L) &&
          (refiner < 16)) {
        refiner++;
        goto justforfun3;
      }
    }

    if (deroundee < fixdigit) {
      correctdigit = correctdigit - (1.0 / pow(10.0, (long double)(j)));
    }
  zerocase3:
    prevtval = prevtval + correctdigit;

    j++;
  }
  t3 = prevtval;
  printf("%s %.15Lf\n", "t3 = ", t3);
  int defaulter;
  long double p[4] = {t0, t1, t2, t3};
  defaulter = 4;
  long double slopeshift = 0;
  int gp = 10;
  printf("\nZero[%.17Lf] = %.15Lf ", t0, Zeroroot(t0, aero));
  printf("\nZero[%.17Lf] = %.15Lf ", t1, Firstroot(t1, aero, xL, yU));
  printf("\nZero[%.17Lf] = %.15Lf ", t2, Secondroot(t2, aero, xL, yB));
  printf("\nZero[%.17Lf] = %.15Lf ", t3, Thirdroot(t3, aero));
  printf("\nt2 / t1 = %.15Lf\n", t2 / t1);
  printf("sqrt (t1/t0) \u2248 %.2Lf\n", sqrtl(fabsl(t1 / t0)));
  printf("sqrt (t2/t0) \u2248 %.2Lf\n", sqrtl(fabsl(t2 / t0)));
  printf("sqrt (t2/t1) \u2248 %.2Lf\n", sqrtl(fabsl(t2 / t1)));

  Tarray *array = malloc(sizeof(Tarray));

  if ((FileProvided && op59 != 14) || (!FileProvided)) {
    printf("Please enter the fineness of the mesh, mp: ");
    scanf("%d", &Mesh1.mp);
    printf("Please enter the numerator of r, r_num:");
    scanf("%d", &Mesh1.r_num);
    printf("Please enter the denomininator or r, r_dem: ");
    scanf("%d", &Mesh1.r_dem);
    printf("Please enter gp (as an integer, perhaps larger than mp): ");
    scanf("%d", &Mesh1.gp);
    printf("Please enter c-parameter, the ratio of the outward seperation "
           "(recommended is 0.123): ");
    scanf("%Lf", &Mesh1.Cpar);
    printf("Please enter trail-parameter, the number of mp's worth of points "
           "til the trailing edge: ");
    scanf("%d", &Mesh1.Tpar);
  };
  r_num = Mesh1.r_num;
  r_dem = Mesh1.r_dem;
  mp = Mesh1.mp;
  gp = Mesh1.gp;
  r = (long double)(((long double)(r_num)) / ((long double)(r_dem)));

  Tvalues(mp, r, p, array, &defaulter, gp);

  long double vi = 24;
  /*
    long double angle;
    printf("Please enter the angle of attack, in degrees less than 90: ");
    scanf("%Lf", &angle);
    long double pi = 3.141592653589;
    FILE *velocity;

    velocity = fopen("U", "w");
    fprintf(velocity, "FoamFile\n");
    fprintf(velocity, "{\n");
    fprintf(velocity, "\tformat      ascii;\n");
    fprintf(velocity, "\tclass       volVectorField;\n");
    fprintf(velocity, "\tobject      U;\n");
    fprintf(velocity, "}\n");
    fprintf(velocity, "// * * * * * * * * * * * * * * * * * * * * * * * * * * *
    * * * * * //\n\n"); fprintf(velocity, "dimensions      [0 1 -1 0 0 0
    0];\n\n"); fprintf(velocity, "\tinternalField   uniform (%.5Lf %.5Lf
    0);\n\n", 24 * cosl(pi * angle / 180.0), 24 * sinl(pi * angle / 180.0));
    fprintf(velocity, "boundaryField\n");
    fprintf(velocity, "{\n");
    fprintf(velocity, "    inlet\n");
    fprintf(velocity, "    {\n");
    fprintf(velocity, "        type            freestreamVelocity;\n");
    fprintf(velocity, "        freestreamValue $internalField;\n");
    fprintf(velocity, "    }\n\n");
    fprintf(velocity, "    outlet\n");
    fprintf(velocity, "    {\n");
    fprintf(velocity, "        type            freestreamVelocity;\n");
    fprintf(velocity, "        freestreamValue $internalField;\n");
    fprintf(velocity, "    }\n\n");
    fprintf(velocity, "    walls\n");
    fprintf(velocity, "    {\n");
    fprintf(velocity, "        type            noSlip;\n");
    fprintf(velocity, "    }\n\n");
    fprintf(velocity, "    frontAndBack\n");
    fprintf(velocity, "    {\n");
    fprintf(velocity, "        type            empty;\n");
    fprintf(velocity, "    }\n");
    fprintf(velocity, "}\n");
  */

  if (aero.m > 0.02) {
    slopeshift = aero.m * 10;
  };
  foilpoint *foil1 = malloc(sizeof(foilpoint));
  foilpoint *foilzoinks = malloc(sizeof(foilpoint));
  trace *trace1 = malloc(sizeof(trace));

  FILE *points[4];
  FILE *faces[3];
  FILE *boundary[2];
  FILE *cells[3];
  FILE *owner[2];
  FILE *neighbour[2];
  headerfill();
  int ncrit = 0;
  int backtrace;
  long double cpar;

  backtrace = Mesh1.Tpar * mp;
  cpar = Mesh1.Cpar / (long double)Mesh1.mp;

  arcadjust(array, aero, &ncrit, &defaulter);
  foilcompute(points, foil1, aero, array, t3, backtrace, trace1, xR);

  printf("\nMesh Data:\n");
  printf("xL: %.5Lf\t%.5Lf\n", xL, Mesh1.xL);
  printf("xR: %.5Lf\t%.5Lf\n", xR, Mesh1.xR);
  printf("yU: %.5Lf\t%.5Lf\n", yU, Mesh1.yU);
  printf("yB: %.5Lf\t%.5Lf\n", yB, Mesh1.yB);
  printf("cpar: %.15Lf\n", cpar);
  printf("m: %.5Lf\t\n", aero.m);
  printf("p: %.5Lf\t\n", aero.p);
  printf("T: %.5Lf\t\n", aero.T);
  printf("t0: %.15Lf\n", t0);
  printf("t1: %.15Lf\n", t1);
  printf("t2: %.15Lf\n", t2);
  printf("t3: %.15Lf\n", t3);
  printf("backtrace: %d\n", backtrace);
  printf("slopeshift: %.5Lf\n", slopeshift);

  printf("ncrit before =%d\n", ncrit);
  OPTIONB opty = four;

  linestack(points, foil1, foilzoinks, aero, cpar, *array, ncrit, Mesh1.xL, mp,
            Mesh1.r_num, Mesh1.r_dem, Mesh1.xR, Mesh1.yU, Mesh1.yB, backtrace,
            trace1, slopeshift, faces, cells, owner, neighbour, boundary, opty);
}
