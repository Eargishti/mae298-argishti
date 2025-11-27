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
  int M0[MAX_SIZE];
  SaveFileMatrixData(matrixfiles[0], matrices, MatrixID, argv[1], stringmats,
                     &count, M0);
  ID = count;
  fclose(matrixfiles[0]);
  matrixfiles[0] = fopen(argv[1], "a");
  Matrix Prod;
  Matrix V_inv;
  Matrix Final;
  Matrix Tran1;
  int nele = 4;
  int nnodes = 4;
  // Properties Matrix
  // A	Izz	  Iyy	J	E	nu	Beta

  // In a 2 dimensional system, the only free DOF at a pin is rotation
  // in fixity, 1 is for free, 0 is for fixed
};
