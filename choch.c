#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void checkface(int ID, FILE *faces, FILE *points, FILE *output) {
  faces = fopen("facecheck", "r");
  points = fopen("figarosafe.txt", "r");
  output = fopen("resultant.txt", "w");

  if (points == NULL || faces == NULL) {
    printf("failed to open file\n");
  }
  int dum;
  int p[4];
  int flabel;
  int bc;
  long double x[4];
  long double y[4];
  long double z[4];

  for (int i = 0; i <= ID; i++) {
    fscanf(faces, "%d %d %d %d %d %d", &p[0], &p[1], &p[2], &p[3], &bc, &flabel);
  };
  for (int k = 0; k < 4; k++) {

    for (int i = 0; i <= p[k]; i++) {
      fscanf(points, "%Lf %Lf %Lf %d", &x[k], &y[k], &z[k], &dum);
    };
    printf("point: %d:  %.15Lf\t%.15Lf\t%.15Lf\n", dum, x[k], y[k], z[k]);
    fprintf(output, "%Lf\t%Lf\n", x[k], y[k]);
    fclose(points);
    points = fopen("figarosafe.txt", "r");
  };

  printf("Average\n (%.15Lf, %.15Lf)\n", 0.25 * (x[0] + x[1] + x[2] + x[3]), 0.25 * (y[0] + y[1] + y[2] + y[3]));
};

void checkall(FILE *faces, FILE *points, FILE *output, int amount, int startFace) {
  faces = fopen("facecheck", "r");
  points = fopen("figarosafe.txt", "r");
  output = fopen("resultant.txt", "w");

  if (points == NULL || faces == NULL || output == NULL) {
    printf("failed to open file\n");
  }
  int dum;
  int p[4];
  int flabel;
  int bc;
  long double x[4];
  long double y[4];
  long double z[4];

  int ID = startFace;
  while (ID < amount + startFace) {

    for (int i = 0; i <= ID; i++) {

      fscanf(faces, "%d %d %d %d %d %d", &p[0], &p[1], &p[2], &p[3], &bc, &flabel);
    };
    for (int k = 0; k < 4; k++) {

      for (int i = 0; i <= p[k]; i++) {
        fscanf(points, "%Lf %Lf %Lf %d", &x[k], &y[k], &z[k], &dum);
      };
      // fprintf(output, "%.15Lf\t%.15Lf\n", 0.25 * (x[0] + x[1] + x[2] + x[3]), 0.25 * (y[0] + y[1] + y[2] + y[3]));
      fclose(points);
      points = fopen("figarosafe.txt", "r");
    };
    fprintf(output, "%.15Lf\t%.15Lf\n", 0.25 * (x[0] + x[1] + x[2] + x[3]), 0.25 * (y[0] + y[1] + y[2] + y[3]));
    fclose(faces);
    faces = fopen("facecheck", "r");

    ID++;
  };
};

int main() {
  int yes = 1;
  int amount = 1;
  printf("Check all( 1 or 0) ?");
  scanf("%d", &yes);
  if (yes) {
    printf("\nnFaces : ");
    scanf("%d", &amount);
  };
  printf("\nPlease enter the label of the face you would like to check: ");
  int label1;
  scanf("%d", &label1);
  FILE *test[3];
  if (yes == 1) {

    checkall(test[0], test[1], test[2], amount, label1);
  } else
    checkface(label1, test[0], test[1], test[2]);
};
