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
  Matrix V_inv;
  Matrix Final;
  Matrix Tran1;
  
	double dx = 1/((double)n-1);
	double x[n] = {0.0};
	Matrix upp = {"upp", {0.0}, n, 1};
	Matrix u = {"solution", {0.0}, n,n};
	StringMatrix sol;

	sol.rows = n;
	sol.columns = 1;

	upp.Element[0][0] = 0;
	u.Element[0][0] = 1.0;

	for (int i = 1; i < (int) n - 1; i++){

	x[i] = ((double) i)*dx;
	upp.Element[i][0] = -sin(M_PI*x[i]) * dx * dx;
	u.Element[i][i - 1] = 1;
	u.Element[i][i] = -2;
	u.Element[i][i + 1] = 1;
	};

	x[(int) n -1] = 1;
	upp.Element[(int)n - 1][0] = 0.0;
	u.Element[(int)n - 1][(int) n - 1] = 1.0;

	Matrix C = MatrixVectorSolve(&u, &sol, &upp);
	FILE *points = fopen("points.txt", "w");
	for (int i = 0; i < (int) n; i++){
	fprintf(points, "%.6lf\t%.6lf\n", x[i], C.Element[i][0]);
	};
	
	



};
