#include "matrix.h"

#define n 10

int main(int argc, char *argv[]) {
  //  infoprint();
  int ID;
  int *MatrixID;
  MatrixID = &ID;
  *MatrixID = 0;
  Matrix *matrices;
  StringMatrix *stringmats;
  MatrixType m11 = Numerical;

  Matrix *temp;
  temp = malloc(sizeof(Matrix));

  Matrix U;

  matrices = malloc(MAX_NUMBER_OF_MATRICES * sizeof(Matrix));

  stringmats = malloc(MAX_NUMBER_OF_MATRICES * sizeof(StringMatrix));

  FILE **matrixfiles;
  switch (argc) {
  case 1:
    printf("Please include an input file\n");
    printf("Usage: ./faytrix [file1] [file2] ...\n");
    return 1;
  default:
    break;
  };
  AllocateFiles(argc, argv, &matrixfiles);
  int count = 0;
  int M0[MAX_SIZE];
  SaveFileMatrixData(matrixfiles[0], matrices, MatrixID, argv[1], stringmats,
                     &count, M0);
  ID = count;
  fclose(matrixfiles[0]);
  matrixfiles[0] = fopen(argv[1], "a");
  FILE *header;
  Matrix Prod;
  Matrix Final;
  Matrix Tran1;
  int N = 2;
  if (argc >= 3)
	  N = atoi(argv[2]);

  /*Matrix M = {"M", {0.0}, N, N};
  
  for (int i = 0; i < N; i++){
  M.Element[i][i] = 1.0/3.0;
     for (int j = i + 1; j < N; j++){
        M.Element[i][j] = 1.0 * (j - i) /8.0;
	 };
  };
 SymmetricUpperT(&M);
 V_inv = inverse(&M); */
long double sum = 0;
FILE *grap = fopen("graphs.txt", "w");

for (int oo = 0; oo < 196; oo++){
  

  LongMatrix M = {"M", {0.0}, N, N};
  
  for (int i = 0; i < N; i++){
  M.Element[i][i] = 1.0;
     for (int j = i + 1; j < N; j++){
        M.Element[i][j] = 0.375 / ((long double)j - (long double) i);
	 };
  };
  LongMatrix V_inv;
 SymmetricUpperLongT(&M);
 V_inv = inverseLong(&M);



  for (int i = 0; i < N; i++){
	  for (int j = 0; j < N; j++){
        sum += V_inv.Element[i][j];
	  };
  };
fprintf(grap, "%.6Lf\t%.8Lf\n", N/sum, logl((long double) N));
N++;
sum = 0;
};


};
