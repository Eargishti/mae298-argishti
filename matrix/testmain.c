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
  FILE *header;
  Matrix Prod;
  Matrix V_inv;
  Matrix Final;
  Matrix Tran1;
  int nele = 4;
  int nnodes = 4;
  // Properties Matrix
  // A	Izz	  Iyy	J	E	nu	Beta
  // In a 2 dimensional system, the only free DOF at a pin is rotation

  // MatrixEnumForm(header, matrices, stringmats, MatrixID, M0);

  // MatrixDefineForm(header, matrices, stringmats, &ID, M0);
#include "defsMatrix.h"
  nnodes = matrices[coord_info].rows;
  nele = matrices[properties].rows;
  StringMatrix u;
  strcpy(u.name, "{u}");
  u.columns = 1;
  u.rows = 6 * nnodes;
  AssignVariables(&u, nnodes);

  AssembleSystemStiffnessMatrix(&matrices[coord_info], &matrices[fixity],

                                &matrices[properties], &matrices[ends], &u,
                                &matrices[concen]);
  m11 = Symbolic;

  // fMatrixPrint(NULL, &u, &m11, stdout);
};
