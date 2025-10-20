#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
// discord.gg/wRBqM5KCTA

#define MAX_SIZE 50
#define MAX_NUMBER_OF_MATRICES 50
#define BUFFER_SIZE 256
#define tol 1E-15
#define VARIABLE_SIZE 6

typedef enum { variable_name, matrix_type, matrix_values } Reading_Type;

typedef enum { Brackets, NoBrackets } mode1;

typedef enum { Numerical, Symbolic } MatrixType;

typedef struct {

  char name[BUFFER_SIZE];
  double Element[MAX_SIZE][MAX_SIZE];
  uint8_t rows;
  uint8_t columns;

} Matrix;

typedef struct {
  char c[VARIABLE_SIZE + 1];
} char3;

typedef struct {
  char name[BUFFER_SIZE];
  char3 Variables[MAX_SIZE][MAX_SIZE];

} StringMatrix;

union Elements {
  double d;
  char3 c;
};

void infoprint() {
  printf("\nWelcome: Program specifications:\n");
  printf("Max n x n of each matrix: %d x %d\n", MAX_SIZE, MAX_SIZE);
  printf("Maximum number of stored matrices: %d\n", MAX_NUMBER_OF_MATRICES);
  printf("Maximum length of each matrix name: %d\n", BUFFER_SIZE - 1);
  printf("Sizes can be changed by modifying the #define statements in the "
         "source file\n");
  printf("\nAbsolute minimum recommended ram size: %lu MB\n",
         10 * MAX_NUMBER_OF_MATRICES * (sizeof(Matrix)) / 1000000);
  printf("float size: %lu, double size: %lu, long double size: %lu\n",
         sizeof(float), sizeof(double), sizeof(long double));
};

void AllocateFiles(int argc, char *argv[], FILE ***matrixfiles) {
  if (argc <= 1) {
    printf("\nAllocateFiles read the file count as %d, and can't allocate %d "
           "files\n",
           argc - 1, argc - 1);
    return;
  };

  *matrixfiles = malloc((argc - 1) * sizeof(FILE *));
  if (*matrixfiles == NULL) {

    printf("malloc failed \n");
    exit(1);
  };

  for (int i = 1; i < argc; i++) {

    (*matrixfiles)[i - 1] = fopen(argv[i], "r");

    if (!(*matrixfiles)[i - 1]) {
      perror(argv[i]);
    };
  };
};

void PrintMatrixData(Matrix *matrices, int k) {

  printf("Variable name: |%s|\n", matrices[k].name);
  printf("%d x %d\n", matrices[k].rows, matrices[k].columns);
  for (int i = 0; i < matrices[k].rows; i++) {
    for (int j = 0; j < matrices[k].columns; j++) {
      printf("%lf\t", matrices[k].Element[i][j]);
    };
    printf("\n");
  };
  printf("\n");
};

void PRintMatrixData(Matrix *matrices, int k) {
  printf("Variable name: |%s|\n", matrices[k].name);
  printf("%d x %d\n", matrices[k].rows, matrices[k].columns);
  for (int i = 0; i < matrices[k].rows; i++) {
    printf("|\t");
    for (int j = 0; j < matrices[k].columns; j++) {
      printf("%lf\t", matrices[k].Element[i][j]);
    };
    printf("|");
    printf("\n");
  };
  printf("\n");
};

void PrintFileContent(FILE **matrixfiles, int filecount) {

  for (int i = 0; i < filecount; i++) {
    FILE *f = matrixfiles[i];
    if (f == NULL) {
      printf("File %d is NULL\n", i);
      continue;
    }

    int ch;
    while ((ch = fgetc(f)) != EOF) {
      putchar(ch);
    }

    printf("\n--- End of file %d ---\n", i);
    rewind(f); // optional, if you want to reuse the file later
  }
}

int IsLetter(char a) {
  if ((a >= 65 && a <= 90) || (a >= 92 && a <= 122) || (a == '_')) {
    return 1;
  } else
    return 0;
};

int IsNumber(char a) {
  if (a >= 48 && a <= 57) {
    return 1;

  } else
    return 0;
  if (a == '\\') {
    printf("Read the %c character\n", a);
    return 0;
  };
};

int IsWhiteSpace(char a) {
  switch (a) {
  case ' ':
  case '\t':
  case '\v':
  case '\f':
  case '\r':
    return 1;
  default:
    return 0;
  };
}

void GetRowData(char *Buffer, Matrix *matrix1, FILE **matrixfile, int row,
                mode1 mode, MatrixType *m1) {
  int BracketFound[2] = {0};
  int i = 0;
  int valuesread = 0;
  int whitespaces = 1;
  char *endptr;
  char *cmprptr;
  if (mode == Brackets && *m1 == Numerical) {
  Reread:
    whitespaces = 1;
    valuesread = 0;
    for (i = 0; i < strlen(Buffer); i++) {
      if (IsWhiteSpace(Buffer[i])) {
        whitespaces++;
        continue;
      };

      if (Buffer[i] == '[' || Buffer[i] == '|') {

        BracketFound[0] = 1;

        for (int j = 0; j < matrix1->columns; j++) {
          cmprptr = j == 0 ? &Buffer[i + 1] : endptr;
          matrix1->Element[row][j] =
              strtod(j == 0 ? &Buffer[i + 1] : endptr, &endptr);
          if (cmprptr != endptr)
            valuesread++;
        };
        if (valuesread == matrix1->columns)
          return;
      };
    };

    if (whitespaces == strlen(Buffer) || valuesread == 0) {
      fgets(Buffer, BUFFER_SIZE, *matrixfile);
      goto Reread;
    };
    if (valuesread != matrix1->columns) {
      printf("In %s, failed to read %d values\n", Buffer, matrix1->columns);
      exit(1);
    };
  }
  if (mode == NoBrackets && *m1 == Numerical) {
  Rer:
    whitespaces = 1;
    valuesread = 0;
    for (i = 0; i < strlen(Buffer); i++) {

      BracketFound[0] = 1;

      for (int j = 0; j < matrix1->columns; j++) {
        cmprptr = j == 0 ? &Buffer[i + 1] : endptr;
        matrix1->Element[row][j] =
            strtod(j == 0 ? &Buffer[i + 1] : endptr, &endptr);
        if (cmprptr != endptr)
          valuesread++;
      };
      if (valuesread == matrix1->columns)
        return;
    };

    if (valuesread == 0) {
      fgets(Buffer, BUFFER_SIZE, *matrixfile);
      goto Rer;
    };
  };
  if (mode == Brackets && *m1 == Symbolic) {
  Rerd:
    whitespaces = 1;
    valuesread = 0;
    for (i = 0; i < strlen(Buffer); i++) {
      if (IsWhiteSpace(Buffer[i])) {
        whitespaces++;
        continue;
      };
    };

    if (valuesread == 0) {
      fgets(Buffer, BUFFER_SIZE, *matrixfile);
      goto Rerd;
    };
  };

  if (valuesread != matrix1->columns) {
    printf("In %s, failed to read %d values\n", Buffer, matrix1->columns);
    exit(1);
  };
}

;

void GetVariableName(char *string1, Matrix *matrix1, char *filename,
                     int *linenumber, FILE *matrixfile, MatrixType *m1) {
  int len = strlen(string1);
  char c;
  int i = 0;
  int j = 0;
  int WasThereEqual = 0;
  int EqualPosition = 0;
  int whitespaces = 0;

  c = string1[0];
  if (c == '$') {
    *m1 = Symbolic;
  } else
    *m1 = Numerical;

  if (!IsLetter(c) && !IsWhiteSpace(c) && c != '\n') {
    printf("c = %d\n", c);
    printf("Isletter(c) = %d, IsWhiteSpace(c) = %d\n", IsLetter(c),
           IsWhiteSpace(c));
    printf("Variable name %s must start with a letter.\n", string1);
    exit(1);
  };

Reread2:
  whitespaces = 0;

  for (i = 0; i < strlen(string1); i++) {
    c = string1[i];

    if (IsWhiteSpace(c) || c == '\n') {
      whitespaces++;
      continue;
    };

    if (!IsLetter(c) && !IsNumber(c) && !IsWhiteSpace(c) || c == '\\') {

    } else {
      matrix1->name[j] = c;
      j++;
    }

    if (string1[i] == '=') {
      WasThereEqual = 1;
      EqualPosition = i;
      matrix1->name[j] = '\0';
      return;
      break;
      goto edge1;
    };
  };

  if (whitespaces == strlen(string1) &&
      fgets(string1, BUFFER_SIZE, matrixfile)) {
    goto Reread2;
  };

edge1:

  if (WasThereEqual == 0) {
    printf("\nIn File %s, line %d, The Variable name %s didn't contain an "
           "equals sign\n",
           filename, *linenumber, string1);
    exit(1);
  };
};

void GetToType(FILE *matrixfile, Matrix *matrices, int MatrixID, char *Buffer) {
  char c[2];
  int i = 0;
  int xWasFound = 0;
  matrices->rows = atoi(Buffer);
  int whitespaces = 0;

  if (matrices->rows <= 0 && Buffer[0] != '\n' && !IsWhiteSpace(Buffer[0])) {
    printf("The type string %s isnt in the format n x m, where n >0\n", Buffer);
    exit(1);
  };

Reread3:
  whitespaces = 0;
  for (i = 0; i < strlen(Buffer); i++) {
    if (Buffer[i] == '\n' || IsWhiteSpace(Buffer[i])) {
      whitespaces++;
    };

    if (Buffer[i] == 'x' || Buffer[i] == 'X') {

      xWasFound = 1;
      matrices->columns = atoi(&Buffer[i + 1]);
      matrices->rows = atoi(&Buffer[0]);
      if (matrices->columns <= 0) {
        printf("The type string %s isn't in the format n x m, where n,m>0\n",
               Buffer);
        exit(1);
      };
      return;
    };
  };
  if (whitespaces == strlen(Buffer) && fgets(Buffer, BUFFER_SIZE, matrixfile)) {
    goto Reread3;
  };

  if (!xWasFound) {
    printf("In format string %s, the character x or X wasnt found\n", Buffer);
    exit(1);
  } else {
  };
};

void SaveFileMatrixData(FILE *matrixfile, Matrix *matrices, int *MatrixID,
                        char *filename) {
  Reading_Type R = variable_name;
  mode1 mode = NoBrackets;
  MatrixType m1 = Numerical;
  printf(mode == NoBrackets ? "File Reading Mode: NoBrackets\n"
                            : "File Reading Mode: With Brackets\n");

  char Buffer[BUFFER_SIZE];
  int ln = 0;
  int *linenumber;
  linenumber = &ln;

  int k = 0;

  // matrices[2]->rows = 4;
  while (fgets(Buffer, BUFFER_SIZE, matrixfile)) {

    GetVariableName(Buffer, &matrices[k], filename, linenumber, matrixfile,
                    &m1);

    fgets(Buffer, BUFFER_SIZE, matrixfile);

    GetToType(matrixfile, &matrices[k], 2, Buffer);

    fgets(Buffer, BUFFER_SIZE, matrixfile);

    for (int i = 0; i < matrices[k].rows; i++) {
      GetRowData(Buffer, &matrices[k], &matrixfile, i, mode, &m1);
      if (i != matrices[k].rows - 1) {
        fgets(Buffer, BUFFER_SIZE, matrixfile);
      };
    };
    printf("Variable name: |%s|\n", matrices[k].name);
    printf("%d x %d\n", matrices[k].rows, matrices[k].columns);
    for (int i = 0; i < matrices[k].rows; i++) {
      for (int j = 0; j < matrices[k].columns; j++) {
        printf("%lf\t", matrices[k].Element[i][j]);
      };
      printf("\n");
    };
    printf("\n");
    k++;
  };
};

double Det(const Matrix *matrix, Matrix *U) {
  // LU Decomposition

  if (matrix->rows != matrix->columns) {
    printf("Determinant not defined for non-square matrices. Matrix %s is %d x "
           "%d\n",
           matrix->name, matrix->rows, matrix->columns);
    return 0.0f;
  };

  int pivot[2] = {0};

  // memcpy(U, matrix, sizeof(Matrix));
  *U = *matrix;
  PrintMatrixData(U, 0);
  double row[MAX_SIZE] = {0};
  double det = 1;
  int k = 0;
  int nonzerofound = 0;
  int pivotwasfound = 0;
  double cof = 1;
  double product[matrix->columns];
  int i = 0;
  int j = 0;

  for (j = 0; j < matrix->columns; j++) {
    pivot[0] = j;

    for (i = j; i < matrix->rows; i++) {
      if (fabs(U->Element[i][j]) <= tol) {
        pivot[0] += !nonzerofound;
        continue;

      } else {
        if (pivot[0] != j && nonzerofound == 0 && pivotwasfound == 0) {
          memcpy(row, U->Element[j], MAX_SIZE * sizeof(double));
          memcpy(U->Element[j], U->Element[pivot[0]],
                 MAX_SIZE * sizeof(double));
          memcpy(U->Element[pivot[0]], row, MAX_SIZE * sizeof(double));
          det *= -1;
          nonzerofound = 1;
          pivotwasfound = 1;
          continue;
        };

        i += (!pivotwasfound);
        pivotwasfound = 1;
        cof = U->Element[i][j] / U->Element[j][j];
        for (k = j; k < matrix->columns; k++) {
          // printf("cof = %.6lf / %.6lf\n", U->Element[i][j],
          // U->Element[j][j]);
          // cof = U->Element[i][k] / U->Element[j][k];
          U->Element[i][k] -= cof * U->Element[j][k];
        };
      };
    }
    product[j] = U->Element[j][j];

    nonzerofound = 0;
    pivotwasfound = 0;
    cof = 1;
  };
  for (int o = 0; o < matrix->columns; o++) {
    det *= product[o];
  };

  printf("\nDet(%s) = %.7lf\n", matrix->name, det);

  return det;
};

double Tran(Matrix matrix, Matrix *T) {
  int i = 0;
  int j = 0;
  *T = matrix;
  T->rows = matrix.columns;
  T->columns = matrix.rows;

  for (i = 0; i < matrix.columns; i++) {
    for (j = 0; j < matrix.rows; j++) {
      T->Element[i][j] = matrix.Element[j][i];
    };
  };
  printf("\n\n");
  PrintMatrixData(T, 0);

  return 0.0f;
};

double Inverse(Matrix *matrix, Matrix *Inv) { return 0.0f; };

void Add(Matrix *A, Matrix *B) {
  Matrix C;
  if (A->rows != B->rows || A->columns != B->columns) {
    printf("\n\nMatrix addition is not defined for matrices of different "
           "size\nMatrix A: %d x %d\nMatrix B: %d x %d\n\n",
           A->rows, A->columns, B->rows, B->columns);
    return;
  };
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      C.Element[i][j] = A->Element[i][j] + B->Element[i][j];
    };
  };
  C.name[0] = '\0';
  printf("\n\n%s + %s = \n", A->name, B->name);
  PrintMatrixData(&C, 0);
};

void Subtract(Matrix *A, Matrix *B) {

  Matrix C;
  if (A->rows != B->rows || A->columns != B->columns) {
    printf("\n\nMatrix Subtraction is not defined for matrices of different "
           "size\nMatrix A: %d x %d\nMatrix B: %d x %d\n\n",
           A->rows, A->columns, B->rows, B->columns);
    return;
  };
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      C.Element[i][j] = A->Element[i][j] - B->Element[i][j];
    };
  };
  C.name[0] = '\0';
  printf("\n\n%s - %s = \n", A->name, B->name);
  PrintMatrixData(&C, 0);
};

void Multiply(Matrix *A, Matrix *B) {
  Matrix C;
  if (A->columns != B->rows) {
    printf("Columns of Matrix A not Equal to rows of Matrix B, can't multiply "
           "(%d x %d) x (%d x %d)\n",
           A->rows, A->columns, B->rows, B->columns);

    return;
  };
  int i = 0;
  int j = 0;
  C.rows = A->rows;
  C.columns = B->columns;
  int i1 = 0;
  int j1 = 0;
  C.name[0] = ' ';
  C.name[1] = '\n';
  C.name[2] = '\0';

  for (j1 = 0; j1 < C.columns; j1++) {

    for (i1 = 0; i1 < C.rows; i1++) {

      for (i = 0; i < B->rows; i++) {

        C.Element[i1][j1] += A->Element[i1][i] * B->Element[i][j1];
      };
    };
  };
  printf("\n\n%s x %s = \n", A->name, B->name);
  PrintMatrixData(&C, 0);

  return;
};

void MatrixVectorSolve() {};

int main(int argc, char *argv[]) {
  //  infoprint();
  int ID;
  int *MatrixID;
  MatrixID = &ID;
  *MatrixID = 0;
  Matrix *matrices;

  Matrix U;
  matrices = malloc(MAX_NUMBER_OF_MATRICES * sizeof(Matrix));

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

  SaveFileMatrixData(matrixfiles[0], matrices, MatrixID, argv[1]);

  Multiply(&matrices[0], &matrices[1]);

  printf("\n\nSize of Matrix: %zu\n\n", sizeof(Matrix));
  printf("\n\nSize of StringMatrix: %zu\n\n", sizeof(StringMatrix));
  /*klocal(1,1) = EAL;
  klocal(1,4) = -EAL;
  klocal(2,2) = k12;
  klocal(2,3) = k6;
  klocal(2,5) = -k12;
  klocal(2,6) = k6;
  klocal(3,3) = k4;
  klocal(3,5) = -k6;
  klocal(3,6) = k2;
  klocal(4,4) = EAL;
  klocal(5,5) = k12;
  klocal(5,6) = -k6;
  klocal(6,6) = k4; */
};
