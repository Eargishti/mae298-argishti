#include "matrix.h"

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
  SaveFileMatrixData(matrixfiles[0], matrices, MatrixID, argv[1], stringmats,
                     &count);

  fclose(matrixfiles[0]);
  matrixfiles[0] = fopen(argv[1], "a");

  Matrix Inv;
  Inverse(&matrices[2], &Inv);
  MatrixVectorSolve(&matrices[2], &stringmats[1], &matrices[3]);
  // SimpleDet(&M1, &U, 1);
};
