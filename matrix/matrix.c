#include "matrix.h"
#include <stdio.h>
#define DIGITS 10

// discord.gg/wRBqM5KCTA

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
// void AssignDOFs(Matrix *matrices

void fMatrixPrint(Matrix *matrices, StringMatrix *stringmatrix, MatrixType *m1,
                  FILE *matrixfile) {

  if (matrixfile == NULL) {
    printf("\nmatrixfile was nulll\n");
    exit(1);
  };

  fprintf(matrixfile, "\n\n");
  if (*m1 == Numerical) {
    fprintf(matrixfile, "%s = \n\n", matrices->name);
    fprintf(matrixfile, "%d x %d\n\n", matrices->rows, matrices->columns);
    fprintf(matrixfile, "\n");
    for (int i = 0; i < matrices->rows; i++) {
      for (int j = 0; j < matrices->columns; j++) {
        fprintf(matrixfile, "%14.7lf\t", matrices->Element[i][j]);
      };
      fprintf(matrixfile, "\n");
    };
  };

  if (*m1 == Symbolic) {
    fprintf(matrixfile, "%s\n", stringmatrix->name);
    fprintf(matrixfile, "%d x %d\n", stringmatrix->rows, stringmatrix->columns);
    fprintf(matrixfile, "\n");
    for (int i = 0; i < stringmatrix->rows; i++) {
      for (int j = 0; j < stringmatrix->columns; j++) {

        fprintf(matrixfile, "%s\t", stringmatrix->Variables[i][j].c);
      };
      fprintf(matrixfile, "\n");
    };
  };
};

void SidebySide(Matrix *matrix1, Matrix *matrix2, FILE *f1) {
  // fprintf(f1, "\n");
  for (int i = 0; i < matrix1->rows; i++) {
    fprintf(f1, "| ");
    for (int j = 0; j < matrix1->columns + matrix2->columns; j++) {
      if (j < matrix1->columns) {
        fprintf(f1, matrix1->Element[i][j] < 0.0 ? "%.4lf  " : "%.5lf ",
                matrix1->Element[i][j]);
      };
      if (j == matrix1->columns) {
        fprintf(f1,
                matrix2->Element[i][j - matrix1->columns] < 0.0 ? "| %.4lf  "
                                                                : "| %.5lf ",
                matrix2->Element[i][j - matrix1->columns]);
      };
      if (j > matrix1->columns) {
        fprintf(f1,
                matrix2->Element[i][j - matrix1->columns] < 0.0 ? "%.4lf  "
                                                                : "%.5lf ",
                matrix2->Element[i][j - matrix1->columns]);
      };
    };
    fprintf(f1, "|\n");
  };
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
                mode1 mode, MatrixType *m1, StringMatrix *stringmat) {
  int BracketFound[2] = {0};
  int i = 0;
  int valuesread = 0;
  int whitespaces = 1;
  int letterfound = 0;
  int sizecounter = 0;
  int k;
  int colindex = 0;

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
        cmprptr = j == 0 ? &Buffer[i] : endptr;
        matrix1->Element[row][j] =
            strtod(j == 0 ? &Buffer[i] : endptr, &endptr);
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

  if (mode == NoBrackets && *m1 == Symbolic) {
  Rerdniga:
    whitespaces = 0;
    valuesread = 0;
    sizecounter = 1;
    colindex = 0;

    for (i = 0; i < strlen(Buffer); i++) {
      k = 0;

      if (IsWhiteSpace(Buffer[i]) || Buffer[i] == '\n') {
        whitespaces++;
        continue;
      };

      if (IsLetter(Buffer[i])) {
        stringmat->Variables[row][colindex].c[k] = Buffer[i];
        k++;
        i++;
        while (IsLetter(Buffer[i]) || IsNumber(Buffer[i]) || Buffer[i] == '_') {
          stringmat->Variables[row][colindex].c[k] = Buffer[i];
          k++;
          i++;
        };
        stringmat->Variables[row][colindex].c[k] = '\0';
        colindex++;
        valuesread++;
        continue;
      };

      if (!IsLetter(Buffer[i])) {
        printf("In string %s, the first character %c wasn't a letter", Buffer,
               Buffer[i]);
        exit(1);
      };
    };

    if (valuesread < stringmat->columns) {
      fgets(Buffer, BUFFER_SIZE, *matrixfile);
      goto Rerdniga;
    };
  };

  if (valuesread != matrix1->columns && valuesread != stringmat->columns) {
    printf("In %s, failed to read %d values\n", Buffer, matrix1->columns);
    exit(1);
  };
}

;

void GetVariableName(char *string1, Matrix *matrix1, char *filename,
                     int *linenumber, FILE *matrixfile, MatrixType *m1,
                     StringMatrix *stringmat1) {
  int len = strlen(string1);
  char c;
  int i = 0;
  int j = 0;
  int WasThereEqual = 0;
  int EqualPosition = 0;
  int whitespaces = 0;
  int firstletterfound = 0;
  char *whichptr;

  c = string1[0];

  if (!IsLetter(c) && !IsWhiteSpace(c) && c != '\n') {
    printf("c = %d\n", c);
    printf("Isletter(c) = %d, IsWhiteSpace(c) = %d\n", IsLetter(c),
           IsWhiteSpace(c));
    printf("Variable name %s must start with a letter.\n", string1);
    exit(1);
  };

Reread2:
  whitespaces = 0;
  firstletterfound = 0;

  for (i = 0; i < strlen(string1); i++) {
    c = string1[i];

    if (IsWhiteSpace(c) || c == '\n') {
      whitespaces++;
      continue;
    };

    if (!IsLetter(c) && !IsNumber(c) && !IsWhiteSpace(c) && c != '$' ||
        c == '\\') {

    } else {
      if (firstletterfound == 0) {
        if (c == '$') {
          *m1 = Symbolic;
        } else {
          *m1 = Numerical;
        };
      };
      //*m1 += (!firstletterfound) * Symbolic * (c == '$');

      firstletterfound = 1;

      matrix1->name[j] = c;
      stringmat1->name[j] = c;
      j++;
    }

    if (string1[i] == '=') {
      WasThereEqual = 1;
      EqualPosition = i;
      matrix1->name[j] = '\0';
      stringmat1->name[j] = '\0';
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

void GetToType(FILE *matrixfile, Matrix *matrices, int MatrixID, char *Buffer,
               StringMatrix *stringmat) {
  char c[2];
  int i = 0;
  int xWasFound = 0;
  matrices->rows = atoi(Buffer);
  stringmat->rows = atoi(Buffer);

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
      stringmat->columns = atoi(&Buffer[i + 1]);
      stringmat->rows = atoi(&Buffer[0]);

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
                        char *filename, StringMatrix *stringmats, int *count) {
  Reading_Type R = variable_name;
  mode1 mode = NoBrackets;
  MatrixType m1 = Numerical;

  printf(mode == NoBrackets ? "File Reading Mode: NoBrackets\n"
                            : "File Reading Mode: With Brackets\n");

  char Buffer[BUFFER_SIZE];
  int ln = 0;
  int *linenumber;
  linenumber = &ln;
  (*count) = 0;
  int k = 0;

  // matrices[2]->rows = 4;
  while (fgets(Buffer, BUFFER_SIZE, matrixfile)) {

    GetVariableName(Buffer, &matrices[k], filename, linenumber, matrixfile, &m1,
                    &stringmats[k]);

    fgets(Buffer, BUFFER_SIZE, matrixfile);

    GetToType(matrixfile, &matrices[k], 2, Buffer, &stringmats[k]);

    fgets(Buffer, BUFFER_SIZE, matrixfile);

    for (int i = 0; i < matrices[k].rows; i++) {
      GetRowData(Buffer, &matrices[k], &matrixfile, i, mode, &m1,
                 &stringmats[k]);
      if (i != matrices[k].rows - 1) {
        fgets(Buffer, BUFFER_SIZE, matrixfile);
      };
    };
    fMatrixPrint(&matrices[k], &stringmats[k], &m1, stdout);
    k++;
    m1 = 0;
    count++;
  };
};

void Multiply(Matrix *A, Matrix *B, Matrix *D) {
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
  *D = C;
  printf("\n\n%s x %s = \n", A->name, B->name);

  fMatrixPrint(&C, NULL, (MatrixType[]){0}, stdout);

  return;
};

double Det(const Matrix *matrix, Matrix *U, int printflag, Matrix *colvector,
           int Swapper[MAX_SIZE][2], int *swapcount) {
  const MatrixType m11 = Numerical;
  // LU Decomposition

  if (matrix->rows != matrix->columns) {
    printf("Determinant not defined for non-square matrices. Matrix %s is %d x "
           "%d\n",
           matrix->name, matrix->rows, matrix->columns);
    return 0.0f;
  } else if (printflag)
    printf("\nRows == columns\n");

  if (colvector->rows != matrix->rows) {
    printf("Column vector in Det function must equivalent number of rows. "
           "%s.rows != %s.rows",
           colvector->name, matrix->name);
    exit(1);
  };
  int pivot[2] = {0};

  // memcpy(U, matrix, sizeof(Matrix));
  *U = *matrix;

  // PrintMatrixData(U, NULL, &m11);
  double row[MAX_SIZE] = {0};
  double det = 1;
  int k = 0;
  int nonzerofound = 0;
  int pivotwasfound = 0;
  double cof = 1;
  double product[matrix->columns];
  int i = 0;
  int j = 0;
  double columnvalue = 1;

  for (j = 0; j < matrix->columns; j++) {
    pivot[0] = j;

    for (i = j; i < matrix->rows; i++) {

      if (fabs(U->Element[i][j]) <= tol) {
        pivot[0] += !nonzerofound;
        continue;

      } else {
        if (pivot[0] != j && nonzerofound == 0 && pivotwasfound == 0) {
          Swapper[*swapcount][0] = j;
          Swapper[*swapcount][1] = pivot[0];
          (*swapcount)++;
          columnvalue = colvector->Element[j][0];
          colvector->Element[j][0] = colvector->Element[pivot[0]][0];
          colvector->Element[pivot[0]][0] = columnvalue;
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
        colvector->Element[i][0] -= cof * colvector->Element[j][0];
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
  /*  printf("Upper Triangular matrix:\n");
    PrintMatrixData(U, NULL, (MatrixType[]){0});
    printf("colvec:\n");
    PrintMatrixData(colvector, NULL, (MatrixType[]){0}); */
  if (printflag)
    printf("\nDet(%s) = %.7lf\n", matrix->name, det);

  return det;
};

double SimpleDet(const Matrix *matrix, Matrix *U, int printflag) {
  const MatrixType m11 = Numerical;
  // LU Decomposition

  if (matrix->rows != matrix->columns) {
    printf("Determinant not defined for non-square matrices. Matrix %s is %d x "
           "%d\n",
           matrix->name, matrix->rows, matrix->columns);
    return 0.0f;
  } else if (printflag)
    printf("\nRows == columns\n");

  int pivot[2] = {0};

  // memcpy(U, matrix, sizeof(Matrix));
  *U = *matrix;

  // PrintMatrixData(U, NULL, &m11);
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
  // printf("Upper Triangular matrix:\n");
  // PrintMatrixData(U, NULL, (MatrixType[]){0});
  if (printflag)
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
  fMatrixPrint(T, NULL, 0, stdout);

  return 0.0f;
};

double Inverse(Matrix *matrix, Matrix *Inv) {
  if (matrix->rows != matrix->columns) {
    printf("Inverse matrix rows aren't columns\n");
  };
  int i = 0;
  int j = 0;
  int k = 0;
  int l = 0;
  Inv->rows = matrix->rows;
  Inv->columns = matrix->columns;

  for (i = 0; i < Inv->rows; i++) {

    Inv->Element[i][i] = 1;
  };
  memcpy(Inv->name, matrix->name, strlen(matrix->name));
  i = strlen(matrix->name);

  Inv->name[i] = '_';
  Inv->name[i + 1] = 'I';
  Inv->name[i + 2] = 'n';
  Inv->name[i + 3] = 'v';
  Inv->name[i + 4] = '\0';
  int pivotwasfound = 0;
  int nonzerofound = 0;
  int pivot = 0;
  double cof = 1;
  Matrix U;
  double det1 = SimpleDet(matrix, &U, 0);
  if (fabs(det1) <= tol) {
    printf("Inverse not defined for 0-determinant matrix\n");
    printf("det(%s) = %.15lf\n", matrix->name, det1);
    exit(1);
  };
  U = *matrix;

  double rowMat[MAX_SIZE] = {0};
  double rowInv[MAX_SIZE] = {0};
  printf("Inverse Matrix before\n");
  for (j = 0; j < matrix->columns; j++) {
    pivot = j;
    pivotwasfound = 0;
    for (int i = j; i < matrix->rows; i++) {
      if (fabs(U.Element[i][j]) <= tol && !pivotwasfound) {
        pivot++;
        continue;
      };
      if (fabs(U.Element[i][j]) >= tol && !pivotwasfound) {
        pivotwasfound = 1;
        if (pivot > j) {
          memcpy(rowMat, U.Element[j], MAX_SIZE * sizeof(double));
          memcpy(rowInv, Inv->Element[j], MAX_SIZE * sizeof(double));
          memcpy(U.Element[j], U.Element[pivot], MAX_SIZE * sizeof(double));
          memcpy(Inv->Element[j], Inv->Element[pivot],
                 MAX_SIZE * sizeof(double));
          memcpy(U.Element[pivot], rowMat, MAX_SIZE * sizeof(double));
          memcpy(Inv->Element[pivot], rowInv, MAX_SIZE * sizeof(double));
        };
        for (l = 0; l < matrix->rows; l++) {
          if (l == j) {
            continue;
          };
          cof = U.Element[l][j] / U.Element[j][j];
          for (k = j; k < U.columns; k++) {
            U.Element[l][k] -= cof * U.Element[j][k];
          };
          for (k = 0; k < Inv->columns; k++) {
            Inv->Element[l][k] -= cof * Inv->Element[j][k];
          };
        };
        goto Invexit;
      };
    };

  Invexit:
  };

  for (i = 0; i < U.rows; i++) {
    cof = U.Element[i][i];
    U.Element[i][i] = 1.0;
    for (j = 0; j < U.columns; j++) {
      Inv->Element[i][j] /= cof;
    };
  };
  fMatrixPrint(matrix, NULL, (MatrixType[]){Numerical}, stdout);
  fMatrixPrint(Inv, NULL, (MatrixType[]){Numerical}, stdout);
  Matrix C;
  Multiply(matrix, Inv, &C);

  return 0.0f;
};

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
  fMatrixPrint(&C, NULL, 0, stdout);
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
  fMatrixPrint(&C, NULL, 0, stdout);
};

void CreateRandomMatrix(Matrix *matrix, unsigned int rows, unsigned int columns,
                        int max, int scale, int *ID) {

  strcpy(matrix->name, "CreatedMatrix");
  if (rows > MAX_SIZE || columns > MAX_SIZE) {
    printf("%d > %d, %d > %d", rows, MAX_SIZE, columns, MAX_SIZE);
    exit(1);
  };
  matrix->rows = rows;
  matrix->columns = columns;

  for (int i = 0; i < matrix->rows; i++) {
    for (int j = 0; j < matrix->rows; j++) {
      matrix->Element[i][j] = (double)(rand() % max) / (double)scale;
      int zc = (rand() % max) < .09 * max;
      int s1 = rand() % 2;
      int g1 = s1 - (s1 + !s1) * !s1;
      matrix->Element[i][j] *= g1 * !zc;
    };
  };
};

void InitializeValues(Matrix *matrix) {
  for (int i = 0; i < matrix->rows; i++) {
    for (int j = 0; j < matrix->columns; j++) {
      matrix->Element[i][j] = 0.0L;
    };
  };
};

void MatrixVectorSolve(Matrix *matrix, StringMatrix *vars, Matrix *cols) {
  if (matrix->rows != matrix->columns) {
    printf("\nNo unique solution exists since matrix %s isn't square, its %d x "
           "%d\n",
           matrix->name, matrix->rows, matrix->columns);
    return;
  };

  if (vars->columns != 1 || cols->columns != 1) {
    printf("Matrix %s or %s is %d x %d or %d x %d and isn't a column vector, "
           "perhaps take its transpose first",
           vars->name, cols->name, vars->rows, vars->columns, cols->rows,
           cols->columns);
    return;
  };

  if (!(cols->rows == vars->rows && vars->rows == matrix->rows)) {
    printf("The rows [%s] {%s} = {%s} aren't equivalennt.\n", matrix->name,
           vars->name, cols->name);
    return;
  };
  int i = 1;
  int j = 1;
  int pivotwasfound = 0;
  int k = 0;
  int pivot = 0;
  double cof = 1;

  Matrix U;
  Matrix C;
  U = *matrix;
  C = *cols;
  int Swapper[MAX_SIZE][2];
  int Swapcount = 0;

  Det(matrix, &U, 0, cols, Swapper, &Swapcount);
  printf("\nUpper Triangular Matrix:\n");
  fMatrixPrint(&U, NULL, (MatrixType[]){0}, stdout);
  printf("\nOriginal column Matrix:\n");
  fMatrixPrint(cols, NULL, (MatrixType[]){0}, stdout);
  for (j = 1; j < U.columns; j++) {
    pivot = j;
    for (i = j - 1; i >= 0; i--) {
      cof = U.Element[i][j] / U.Element[j][j];
      cols->Element[i][0] -= cof * cols->Element[j][0];
      for (k = j; k < U.columns; k++) {
        U.Element[i][k] -= cof * U.Element[j][k];
      };
    };
  };
  for (int o0 = 0; o0 < vars->rows; o0++) {
    C.Element[o0][0] = cols->Element[o0][0] / U.Element[o0][o0];
  };
  /*printf("Swapcount = %d\n", Swapcount);
  for (int oi = 0; oi < Swapcount; oi++) {
    printf("Swapping %s and %s\n", vars->Variables[Swapper[oi][0]][0].c,
           vars->Variables[Swapper[oi][1]][0].c);
    cof = C.Element[Swapper[oi][0]][0];
    C.Element[Swapper[oi][0]][0] = C.Element[Swapper[oi][1]][0];
    C.Element[Swapper[oi][1]][0] = cof;
  };*/

  printf("\nU.columns = %d\tU.rows = %d\n", U.columns, U.rows);
  printf("\nDiagonalized Matrix:\n");
  fMatrixPrint(&U, NULL, (MatrixType[]){0}, stdout);
  printf("\nNew Column Matrix:\n");
  fMatrixPrint(cols, NULL, (MatrixType[]){0}, stdout);

  printf("vars.rows = %d\n", vars->rows);
  for (int ou = 0; ou < vars->rows; ou++) {
    printf("%s = %.18lf\n", vars->Variables[ou][0].c, C.Element[ou][0]);
  };
};
