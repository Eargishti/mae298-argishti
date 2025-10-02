#include <stdio.h>
#include <string.h>

void velcro(FILE *cells, int ID, FILE *faces, FILE *points, FILE *output, FILE *output2) {

  FILE *fp = fopen("cellcheck", "rb");

  if (!fp) {
    printf("\nFile failed to open\n");
    return;
  }

  printf("Enter Cell to find:\n");
  size_t n = ID;
  //
  long long off = (long long)(n) * 7 * sizeof(int);

  if (fseek(fp, off, SEEK_SET) != 0) {
    printf("fseek failed to run\n");
    return;
  }

  int face[7];
  if (fread(face, sizeof(int), 7, fp) != 7)
    return;

  printf("You picked %zu: Cell %d: \n", n, face[6]);
  printf("%d\t%d\t%d\t%d\t%d\t%d\n", face[0], face[1], face[2], face[3], face[4], face[5]);
  fclose(fp);

  faces = fopen("facecheck", "rb");
  points = fopen("pointcheck", "rb");
  output = fopen("hoh", "w");
  output2 = fopen("quick", "w");
  FILE *diffout[6];
  char m = 48;
  char *stong = "faCeC";
  diffout[0] = fopen("f0", "w");
  diffout[1] = fopen("f1", "w");
  diffout[2] = fopen("f2", "w");
  diffout[3] = fopen("f3", "w");
  diffout[4] = fopen("f4", "w");
  diffout[5] = fopen("f5", "w");

  // int p[4][6];
  int p[6][4];
  int flabel[6];

  int bc[6];

  long double x[4][6];
  long double y[4][6];
  long double z[4][6];

  const long long FACE_REC_SIZE = 6LL * sizeof(int);
  const long long POINT_REC_SIZE = 3LL * sizeof(long double) + 1LL * sizeof(int);

  for (int k = 0; k < 6; k++) {
    off = (long long)face[k] * FACE_REC_SIZE;

    fseek(faces, off, SEEK_SET);
    fread(p[k], sizeof(int), 4, faces);

    fread(&bc[k], sizeof(int), 1, faces);
    fread(&flabel[k], sizeof(int), 1, faces);
    printf("Face %d: %d\t%d\t%d\t%d\n", face[k], p[k][0], p[k][1], p[k][2], p[k][3]);
  };

  int skip;
  printf("Cell %d: \n", face[6]);
  for (int k = 0; k < 6; k++) {
    printf("Face %d: %d\tBC: %d\n", k, face[k], bc[k]);
    for (int i = 0; i < 4; i++) {

      for (int j = 0; j <= p[k][i]; j++) {
        off = p[k][i] * POINT_REC_SIZE;
        fseek(points, off, SEEK_SET);
        fread(&x[i][k], sizeof(long double), 1, points);

        fread(&y[i][k], sizeof(long double), 1, points);
        fread(&z[i][k], sizeof(long double), 1, points);
        fread(&skip, sizeof(int), 1, points);
        // fscanf(points, "%Lf %Lf %Lf %d", &x[i][k], &y[i][k], &z[i][k], &skip);
      };
      printf("\t%.15Lf\t%.15Lf\t%.4Lf\n", x[i][k], y[i][k], z[i][k]);
      fprintf(output, "%.15Lf\t%.15Lf\n", x[i][k], y[i][k]);
      fprintf(diffout[k], "%.15Lf\t%.15Lf\n", x[i][k], y[i][k]);
    }
    printf("Average (%.15Lf, %.15Lf)\n", 0.25 * (x[0][k] + x[1][k] + x[2][k] + x[3][k]), 0.25 * (y[0][k] + y[1][k] + y[2][k] + y[3][k]));
    fprintf(diffout[k], "%.15Lf\t%.15Lf\n", 0.25 * (x[0][k] + x[1][k] + x[2][k] + x[3][k]), 0.25 * (y[0][k] + y[1][k] + y[2][k] + y[3][k]));

    // fprintf(output, "%.15Lf\t%.15Lf\n", 0.25 * (x[0][k] + x[1][k] + x[2][k] + x[3][k]), 0.25 * (y[0][k] + y[1][k] + y[2][k] + y[3][k]));
    if (bc[k] == 5 || bc[k] == 6) {
      fprintf(output2, "%.15Lf\t%.15Lf\n", 0.25 * (x[0][k] + x[1][k] + x[2][k] + x[3][k]), 0.25 * (y[0][k] + y[1][k] + y[2][k] + y[3][k]));
    };
  };

  printf("cell %d\n", face[6]);
  printf("\nCell output files named \"hoh\" and \"quick\"\n");
};

int main() {

  FILE *cells;
  int ID;
  FILE *faces;
  FILE *points;
  FILE *output;
  FILE *output2;
  printf("\nEnter cell to find: ");
  scanf("%d", &ID);

  velcro(cells, ID, faces, points, output, output2);
}
/*

  FILE *fp = fopen("cellcheck", "rb");

  if (!fp) {
    printf("\nFile failed to open\n");
    return 1;
  }

  printf("Enter Cell to find:\n");
  size_t n = 1;
  scanf("%zu", &n); // 1-based record index
                    //
  long long off = (long long)(n - 1) * 7 * sizeof(int);

  if (fseek(fp, off, SEEK_SET) != 0) {
    printf("fseek failed to run\n");
    return 1;
  }

  int f[7];
  if (fread(f, sizeof(int), 7, fp) != 7)
    return 1;

  printf("You picked %zu: Cell %d: \n", n, f[6]);
  printf("%d\t%d\t%d\t%d\t%d\t%d\n", f[0], f[1], f[2], f[3], f[4], f[5]);
  fclose(fp);

  faces = fopen("facecheck", "r");
  points = fopen("pointcheck", "r");
  output = fopen("hoh", "w");
  output2 = fopen("quick", "w");

  int p[4][6];

  int flabel[6];

  int bc[6];

  long double x[4][6];
  long double y[4][6];
  long double z[4][6];

  for (int k = 0; k < 6; k++) {

    for (int j = 0; j <= face[k]; j++) {

      for (int i = 0; i < 4; i++) {
        fscanf(faces, "%d", &p[i][k]);
      };

      fscanf(faces, "%d %d", &bc[k], &flabel[k]);
    };
    fclose(faces);
    faces = fopen("facecheck", "r");
  };

  int skip;
  printf("Cell %d: \n", dum);
  for (int k = 0; k < 6; k++) {
    printf("Face %d: %d\tBC: %d\n", k, face[k], bc[k]);
    for (int i = 0; i < 4; i++) {

      for (int j = 0; j <= p[i][k]; j++) {

        fscanf(points, "%Lf %Lf %Lf %d", &x[i][k], &y[i][k], &z[i][k], &skip);
      };
      printf("\t%.15Lf\t%.15Lf\t%.4Lf\n", x[i][k], y[i][k], z[i][k]);
      fprintf(output, "%.15Lf\t%.15Lf\n", x[i][k], y[i][k]);
      fclose(points);
      points = fopen("figarosafe.txt", "r");
    }
    printf("Average (%.15Lf, %.15Lf)\n", 0.25 * (x[0][k] + x[1][k] + x[2][k] + x[3][k]), 0.25 * (y[0][k] + y[1][k] + y[2][k] + y[3][k]));
    // fprintf(output, "%.15Lf\t%.15Lf\n", 0.25 * (x[0][k] + x[1][k] + x[2][k] + x[3][k]), 0.25 * (y[0][k] + y[1][k] + y[2][k] + y[3][k]));
    if (bc[k] == 5 || bc[k] == 6) {
      fprintf(output2, "%.15Lf\t%.15Lf\n", 0.25 * (x[0][k] + x[1][k] + x[2][k] + x[3][k]), 0.25 * (y[0][k] + y[1][k] + y[2][k] + y[3][k]));
    };
  };

  printf("cell %d\n", dum);
}; */
