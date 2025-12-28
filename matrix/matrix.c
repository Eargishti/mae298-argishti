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

void fMatrixPrint(const Matrix *matrices, StringMatrix *stringmatrix,
                  MatrixType *m1, FILE *matrixfile) {

  if (matrixfile == NULL) {
    printf("\nmatrixfile was null\n");
    exit(1);
  };

  fprintf(matrixfile, "\n\n");
  if (*m1 == Numerical) {
    fprintf(matrixfile, "%s = \n\n", matrices->name);
    fprintf(matrixfile, "%d x %d\n\n", matrices->rows, matrices->columns);
    fprintf(matrixfile, "\n");
    for (int i = 0; i < matrices->rows; i++) {
      for (int j = 0; j < matrices->columns; j++) {
        fprintf(matrixfile, "%+4.7lf\t", matrices->Element[i][j]);
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

void PRintMatrixData(const Matrix *matrices) {
  printf("%s = \n", matrices->name);
  printf("%d x %d\n", matrices->rows, matrices->columns);

  for (int i = 0; i < matrices->rows; i++) {
    for (int j = 0; j < matrices->columns; j++) {
      fprintf(stdout, "%.7lf\t", matrices->Element[i][j]);
    };
    fprintf(stdout, "\n");
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

void ShishWhiteSpace(char a, int *value) {
  switch (a) {
  case ' ':
  case '\t':
  case '\v':
  case '\f':
  case '\r':
    *value = 1;
  default:
    *value = 0;
  };
};

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
          firstletterfound = 1;
          continue;
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
                        char *filename, StringMatrix *stringmats, int *count,
                        int m[MAX_SIZE]) {
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
    m[k] = m1;

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
    (*count)++;
  };

  m[k] = -1;
};

Matrix Multiply(const Matrix *A, const Matrix *B) {
  Matrix C = {"", {0.0}, A->rows, B->columns};

  if (A->columns != B->rows) {
    printf("Columns of Matrix A not Equal to rows of Matrix B, can't multiply "
           "(%d x %d) x (%d x %d)\n",
           A->rows, A->columns, B->rows, B->columns);

    return C;
  };
  int i = 0;
  int j = 0;
  C.rows = A->rows;
  C.columns = B->columns;
  int i1 = 0;
  int j1 = 0;
  int o = 0;
  const char gog[] = " x ";
  strcpy(C.name, A->name);
  strcat(C.name, gog);
  strcat(C.name, B->name);

  for (j1 = 0; j1 < C.columns; j1++) {

    for (i1 = 0; i1 < C.rows; i1++) {

      for (i = 0; i < B->rows; i++) {

        C.Element[i1][j1] += A->Element[i1][i] * B->Element[i][j1];
      };
    };
  };
  return C;
};

Matrix DebugMultiply(const Matrix *A, const Matrix *B, int problematicrow[2]) {

  Matrix C = {"", {0.0}, A->rows, B->columns};

  if (A->columns != B->rows) {
    printf("Columns of Matrix A not Equal to rows of Matrix B, can't multiply "
           "(%d x %d) x (%d x %d)\n",
           A->rows, A->columns, B->rows, B->columns);

    return C;
  };
  int i = 0;
  int j = 0;
  C.rows = A->rows;
  C.columns = B->columns;
  int i1 = 0;
  int j1 = 0;
  int o = 0;
  int printmode = 0;

  for (j1 = 0; j1 < C.columns; j1++) {

    for (i1 = 0; i1 < C.rows; i1++) {
      if (i1 == problematicrow[0] || i1 == problematicrow[1]) {
        printmode = 1;
      };

      for (i = 0; i < B->rows; i++) {

        C.Element[i1][j1] += A->Element[i1][i] * B->Element[i][j1];
        if (printmode) {
          printf("%s(%d)(%d) += %.7lf * %.7lf\n", "R", i1, j1,
                 A->Element[i1][i], B->Element[i][j1]);
        };
      };
      printmode = 0;
    };
  };
  return C;
};

SmallMatrix SmallMultiply(const SmallMatrix *A, const SmallMatrix *B) {
  SmallMatrix C = {"", {0.0}, A->rows, B->columns};

  if (A->columns != B->rows) {
    printf("Columns of Matrix A not Equal to rows of Matrix B, can't multiply "
           "(%d x %d) x (%d x %d)\n",
           A->rows, A->columns, B->rows, B->columns);

    return C;
  };
  int i = 0;
  int j = 0;
  C.rows = A->rows;
  C.columns = B->columns;
  int i1 = 0;
  int j1 = 0;
  int o = 0;
  const char gog[] = " x ";
  strcpy(C.name, A->name);
  strcat(C.name, gog);
  strcat(C.name, B->name);

  for (j1 = 0; j1 < C.columns; j1++) {

    for (i1 = 0; i1 < C.rows; i1++) {

      for (i = 0; i < B->rows; i++) {

        C.Element[i1][j1] += A->Element[i1][i] * B->Element[i][j1];
      };
    };
  };
  return C;
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
  const char gog[] = "_Tran";
  strcat(T->name, gog);

  for (i = 0; i < matrix.columns; i++) {
    for (j = 0; j < matrix.rows; j++) {
      T->Element[i][j] = matrix.Element[j][i];
    };
  };
  printf("\n\n");
  fMatrixPrint(T, NULL, (MatrixType[]){Numerical}, stdout);

  return 0.0f;
};

Matrix Transpose(const Matrix *matrix) {
  Matrix T;
  T.rows = matrix->columns;
  T.columns = matrix->rows;
  const char gog[] = "_T";
  strcpy(T.name, matrix->name);
  strcat(T.name, gog);
  for (int i = 0; i < matrix->columns; i++) {
    for (int j = 0; j < matrix->rows; j++) {
      T.Element[i][j] = matrix->Element[j][i];
    };
  };
  return T;
};

SmallMatrix SmallTranspose(const SmallMatrix *matrix) {
  SmallMatrix T;
  T.rows = matrix->columns;
  T.columns = matrix->rows;
  const char gog[] = "_T";
  strcpy(T.name, matrix->name);
  strcat(T.name, gog);
  for (int i = 0; i < matrix->columns; i++) {
    for (int j = 0; j < matrix->rows; j++) {
      T.Element[i][j] = matrix->Element[j][i];
    };
  };
  return T;
};

double Inverse(const Matrix *matrix, Matrix *Inv) {
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

  Inv->name[i - 1] = '_';
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

  return 0.0f;
};

Matrix inverse(const Matrix *matrix) {
  Matrix Inv = {"", {0.0}, matrix->rows, matrix->columns};
  if (matrix->rows != matrix->columns) {
    printf("Inverse matrix rows aren't columns\n");
  };
  int i = 0;
  int j = 0;
  int k = 0;
  int l = 0;

  for (i = 0; i < Inv.rows; i++) {

    Inv.Element[i][i] = 1;
  };
  strcpy(Inv.name, matrix->name);
  const char gog[] = "_Inv";
  strcat(Inv.name, gog);
  int pivotwasfound = 0;
  int nonzerofound = 0;
  int pivot = 0;
  double cof = 1;
  Matrix U;
  double det1 = SimpleDet(matrix, &U, 0);
  if (fabs(det1) <= tol) {
    printf("Inverse not defined for 0-determinant matrix\n");
    //    fMatrixPrint(matrix, NULL, &(MatrixType){Numerical}, stdout);
    printf("det(%s) = %.15lf\n", matrix->name, det1);
    exit(1);
  };
  U = *matrix;

  double rowMat[MAX_SIZE] = {0};
  double rowInv[MAX_SIZE] = {0};
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
          memcpy(rowInv, Inv.Element[j], MAX_SIZE * sizeof(double));
          memcpy(U.Element[j], U.Element[pivot], MAX_SIZE * sizeof(double));
          memcpy(Inv.Element[j], Inv.Element[pivot], MAX_SIZE * sizeof(double));
          memcpy(U.Element[pivot], rowMat, MAX_SIZE * sizeof(double));
          memcpy(Inv.Element[pivot], rowInv, MAX_SIZE * sizeof(double));
        };
        for (l = 0; l < matrix->rows; l++) {
          if (l == j) {
            continue;
          };
          cof = U.Element[l][j] / U.Element[j][j];
          for (k = j; k < U.columns; k++) {
            U.Element[l][k] -= cof * U.Element[j][k];
          };
          for (k = 0; k < Inv.columns; k++) {
            Inv.Element[l][k] -= cof * Inv.Element[j][k];
          };
        };
        goto Invex;
      };
    };

  Invex:
  };

  for (i = 0; i < U.rows; i++) {
    cof = U.Element[i][i];
    U.Element[i][i] = 1.0;
    for (j = 0; j < U.columns; j++) {
      Inv.Element[i][j] /= cof;
    };
  };
  /* Matrix H;
   H = Multiply(matrix, &Inv);
   MatrixType m1 = Numerical;
   fMatrixPrint(&H, NULL, &m1, stdout); */

  return Inv;
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

Matrix Subtract(const Matrix *A, const Matrix *B) {

  Matrix C = {"", {0.0}, A->rows, A->columns};
  C.rows = A->rows;
  C.columns = A->columns;
  if (A->rows != B->rows || A->columns != B->columns) {
    printf("\n\nMatrix Subtraction is not defined for matrices of different "
           "size\nMatrix A: %d x %d\nMatrix B: %d x %d\n\n",
           A->rows, A->columns, B->rows, B->columns);
    exit(1);
  };
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      C.Element[i][j] = A->Element[i][j] - B->Element[i][j];
    };
  };
  /*strcat(C.name, A->name);
  const char gog[] = " - ";
  strcat(C.name, gog);
  strcat(C.name, B->name); */

  strcpy(C.name, A->name);
  PRintMatrixData(&C);
  return C;
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

Matrix MatrixVectorSolve(Matrix *matrix, StringMatrix *vars, Matrix *cols) {
  Matrix C;
  if (matrix->rows != matrix->columns) {
    printf("\nNo unique solution exists since matrix %s isn't square, its %d x "
           "%d\n",
           matrix->name, matrix->rows, matrix->columns);
    return C;
  };

  if (vars->columns != 1 || cols->columns != 1) {
    printf("Matrix %s or %s is %d x %d or %d x %d and isn't a column vector, "
           "perhaps take its transpose first",
           vars->name, cols->name, vars->rows, vars->columns, cols->rows,
           cols->columns);
    return C;
  };

  if (!(cols->rows == vars->rows && vars->rows == matrix->rows)) {
    printf("The rows [%s] {%s} = {%s} aren't equivalennt.\n", matrix->name,
           vars->name, cols->name);
    return C;
  };
  int i = 1;
  int j = 1;
  int pivotwasfound = 0;
  int k = 0;
  int pivot = 0;
  double cof = 1;

  Matrix U;
  U = *matrix;
  C = *cols;
  int Swapper[MAX_SIZE][2];
  int Swapcount = 0;

  Det(matrix, &U, 0, cols, Swapper, &Swapcount);
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

  printf("vars.rows = %d\n", vars->rows);
  for (int ou = 0; ou < vars->rows; ou++) {
    printf("%s = %.18lf\n", vars->Variables[ou][0].c, C.Element[ou][0]);
  };
  return C;
};

void MatrixEnumForm(FILE *header, Matrix *matrices, StringMatrix *stringmats,
                    int *ID, int m[MAX_SIZE]) {
  header = fopen("enumsMatrix.h", "w");
  fprintf(header, "typedef enum {  ");
  int k = 0;
  char name[BUFFER_SIZE];
  while (matrices[k].rows != 0 || stringmats[k].rows != 0) {
    if (m[k] == Numerical) {
      fprintf(header, "%s , ", matrices[k].name);
    };
    if (m[k] == Symbolic) {
      fprintf(header, "%s , ", stringmats[k].name);
    };
    k++;
  };
  fprintf(header, "  }");
  fprintf(header, " NamesofMats;\n");
  fprintf(header, "\n");
  fclose(header);
};

void MatrixDefineForm(FILE *header, Matrix *matrices, StringMatrix *stringmats,
                      int *ID, int m[MAX_SIZE]) {
  header = fopen("defsMatrix.h", "w");
  int k = 0;
  char name[BUFFER_SIZE];

  while (matrices[k].rows != 0 || stringmats[k].rows != 0) {
    if (m[k] == Numerical) {
      fprintf(header, "#define %s %d\n", matrices[k].name, k);
    };
    if (m[k] == Symbolic) {
      fprintf(header, "#define %s %d\n", stringmats[k].name, k);
    };
    k++;
  };
};

void SymmetricUpperT(Matrix *matrix) {
  int i = 0;
  int j = 0;
  for (j = 0; j < matrix->rows; j++) {
    for (i = 0; i < j; i++) {
      matrix->Element[j][i] = matrix->Element[i][j];
    };
  };
};

void SymmetricSmallUT(SmallMatrix *matrix) {
  int i = 0;
  int j = 0;
  for (j = 0; j < matrix->rows; j++) {
    for (i = 0; i < j; i++) {
      matrix->Element[j][i] = matrix->Element[i][j];
    };
  };
};

Matrix elk(double A, double Izz, double Iyy, double J, double E, double nu,
           double L) {
  Matrix Elk;
  char Nom[] = "estiff";
  memcpy(Elk.name, Nom, sizeof(Nom));
  Elk.rows = 12;
  Elk.columns = 12;
  InitializeValues(&Elk);
  Elk.Element[0][0] = A * E / L;
  Elk.Element[0][6] = -A * E / L;
  Elk.Element[1][1] = 12 * Izz * E / (L * L * L);
  Elk.Element[1][5] = 6 * Izz * E / (L * L);
  Elk.Element[1][7] = -12 * Izz * E / (L * L * L);
  Elk.Element[1][11] = 6 * Izz * E / (L * L);
  Elk.Element[2][2] = 12 * Iyy * E / (L * L * L);
  Elk.Element[2][4] = -6 * Iyy * E / (L * L);
  Elk.Element[2][8] = -12 * Iyy * E / (L * L * L);
  Elk.Element[2][10] = -6 * Iyy * E / (L * L);
  Elk.Element[3][3] = J * E / (2 * (1 + nu) * L);
  Elk.Element[3][9] = -J * E / (2 * (1 + nu) * L);
  Elk.Element[4][4] = 4 * Iyy * E / L;
  Elk.Element[4][8] = 6 * Iyy * E / (L * L);
  Elk.Element[4][10] = 2 * Iyy * E / L;
  Elk.Element[5][5] = 4 * Izz * E / L;
  Elk.Element[5][7] = -6 * Izz * E / (L * L);
  Elk.Element[5][11] = 2 * Izz * E / L;
  Elk.Element[6][6] = A * E / L;
  Elk.Element[7][7] = 12 * Izz * E / (L * L * L);
  Elk.Element[7][11] = -6 * Izz * E / (L * L);
  Elk.Element[8][8] = 12 * Iyy * E / (L * L * L);
  Elk.Element[8][10] = 6 * Iyy * E / (L * L);
  Elk.Element[9][9] = J * E / (2 * (1 + nu) * L);
  Elk.Element[10][10] = 4 * Iyy * E / L;
  Elk.Element[11][11] = 4 * Izz * E / L;
  SymmetricUpperT(&Elk);
  return Elk;
};

SmallMatrix Smallelk(double A, double Izz, double Iyy, double J, double E,
                     double nu, double L) {
  SmallMatrix Elk = {"esitff", {0.0}, 12, 12};
  Elk.Element[0][0] = A * E / L;
  Elk.Element[0][6] = -A * E / L;
  Elk.Element[1][1] = 12 * Izz * E / (L * L * L);
  Elk.Element[1][5] = 6 * Izz * E / (L * L);
  Elk.Element[1][7] = -12 * Izz * E / (L * L * L);
  Elk.Element[1][11] = 6 * Izz * E / (L * L);
  Elk.Element[2][2] = 12 * Iyy * E / (L * L * L);
  Elk.Element[2][4] = -6 * Iyy * E / (L * L);
  Elk.Element[2][8] = -12 * Iyy * E / (L * L * L);
  Elk.Element[2][10] = -6 * Iyy * E / (L * L);
  Elk.Element[3][3] = J * E / (2 * (1 + nu) * L);
  Elk.Element[3][9] = -J * E / (2 * (1 + nu) * L);
  Elk.Element[4][4] = 4 * Iyy * E / L;
  Elk.Element[4][8] = 6 * Iyy * E / (L * L);
  Elk.Element[4][10] = 2 * Iyy * E / L;
  Elk.Element[5][5] = 4 * Izz * E / L;
  Elk.Element[5][7] = -6 * Izz * E / (L * L);
  Elk.Element[5][11] = 2 * Izz * E / L;
  Elk.Element[6][6] = A * E / L;
  Elk.Element[7][7] = 12 * Izz * E / (L * L * L);
  Elk.Element[7][11] = -6 * Izz * E / (L * L);
  Elk.Element[8][8] = 12 * Iyy * E / (L * L * L);
  Elk.Element[8][10] = 6 * Iyy * E / (L * L);
  Elk.Element[9][9] = J * E / (2 * (1 + nu) * L);
  Elk.Element[10][10] = 4 * Iyy * E / L;
  Elk.Element[11][11] = 4 * Izz * E / L;
  SymmetricSmallUT(&Elk);
  return Elk;
};

SmallMatrix ShearAdjusteElk(double A, double Izz, double Iyy, double J,
                            double E, double nu, double L, double As2,
                            double As3) {
  SmallMatrix Elk = {"esitff", {0.0}, 12, 12};
  double G = E / (2 * (1 + nu));
  double eta_y = E * Izz / (G * As2);
  double kgamma_y = E * Izz / (L * ((L * L) / 12.0 + eta_y));
  double eta_z = E * Iyy / (G * As3);
  double kgamma_z = E * Iyy / (L * ((L * L) / 12.0 + eta_z));
  Elk.Element[0][0] = A * E / L;
  Elk.Element[0][6] = -A * E / L;
  Elk.Element[1][1] = kgamma_y;
  Elk.Element[1][5] = kgamma_y * (L / 2.0);
  Elk.Element[1][7] = -kgamma_y;
  Elk.Element[1][11] = kgamma_y * (L / 2.0);
  Elk.Element[2][2] = kgamma_z;
  Elk.Element[2][4] = -kgamma_z * (L / 2.0);
  Elk.Element[2][8] = -kgamma_z;
  Elk.Element[2][10] = -kgamma_z * (L / 2.0);
  Elk.Element[3][3] = J * E / (2 * (1 + nu) * L);
  Elk.Element[3][9] = -J * E / (2 * (1 + nu) * L);
  Elk.Element[4][4] = kgamma_z * (L * L / 3.0 + eta_z);
  Elk.Element[4][8] = kgamma_z * (L / 2.0);
  Elk.Element[4][10] = kgamma_z * (L * L / 6.0 - eta_z);
  Elk.Element[5][5] = kgamma_y * (L * L / 3.0 + eta_y);
  Elk.Element[5][7] = -kgamma_y * (L / 2.0);
  Elk.Element[5][11] = kgamma_y * (L * L / 6.0 - eta_y);
  Elk.Element[6][6] = A * E / L;
  Elk.Element[7][7] = kgamma_y;
  Elk.Element[7][11] = -kgamma_y * (L / 2.0);
  Elk.Element[8][8] = kgamma_z;
  Elk.Element[8][10] = kgamma_z * (L / 2.0);
  Elk.Element[9][9] = J * E / (2 * (1 + nu) * L);
  Elk.Element[10][10] = kgamma_z * (L * L / 3.0 + eta_z);
  Elk.Element[11][11] = kgamma_y * (L * L / 3.0 + eta_y);
  SymmetricSmallUT(&Elk);
  return Elk;
};

void norm(double xaxis[3]) {
  double n1 =
      sqrt(xaxis[0] * xaxis[0] + xaxis[1] * xaxis[1] + xaxis[2] * xaxis[2]);
  for (int i = 0; i < 3; i++) {
    xaxis[i] /= n1;
  };
};

Matrix GammaMat(double beta, double xaxis[3]) {
  Matrix etran = {"etran", {0.0}, 12, 12};
  // When looking at node i of the member towards node j, Beta is the clockwise
  // angle downards from the original median
  double ztemp[3] = {0.0};
  double ytemp[3] = {0.0};

  double zaxis[3] = {0.0};
  double yaxis[3] = {0.0};
  norm(xaxis);
  int isGimbalLock = (1.0 - fabs(xaxis[1])) <= 1E+06 * tol;

  if (1.0 - fabs(xaxis[1]) <= 1E+06 * tol) {
    printf("\nGimbal lock\n");
    /*       | i 	j 	k  |
   y' =      | 0	0	1  |
             | x0	x1	x2 | */
    ytemp[0] = -xaxis[1];
    ytemp[1] = xaxis[0];
    ytemp[2] = 0;
    norm(ytemp);

    /*    | i 	j 	k  |
   z =    | x0  x1   x2  |
          | -x1	x0	0  | */
    ztemp[0] = -xaxis[2] * xaxis[0];
    ztemp[1] = -xaxis[2] * xaxis[1];
    ztemp[2] = xaxis[0] * xaxis[0] + xaxis[1] * xaxis[1];
    norm(ztemp);
  } else {
    /* | i 	j 	k  |
       | x0	x1	x2 |
       | 0	1	0  | */

    ztemp[0] = -xaxis[2];
    ztemp[1] = 0;
    ztemp[2] = xaxis[0];
    norm(ztemp);
    /* | i 	j 	k  |
       | -x2  0   x0  |
       | x0	x1	x2  | */

    ytemp[0] = -xaxis[0] * xaxis[1];
    ytemp[1] = xaxis[0] * xaxis[0] + xaxis[2] * xaxis[2];
    ytemp[2] = -xaxis[2] * xaxis[1];
    norm(ytemp);
  };

  zaxis[0] = cosl(beta) * ztemp[0] - sinl(beta) * ytemp[0];
  zaxis[1] = cosl(beta) * ztemp[1] - sinl(beta) * ytemp[1];
  zaxis[2] = cosl(beta) * ztemp[2] - sinl(beta) * ytemp[2];
  norm(zaxis);

  yaxis[0] = cosl(beta) * ytemp[0] + sinl(beta) * ztemp[0];
  yaxis[1] = cosl(beta) * ytemp[1] + sinl(beta) * ztemp[1];
  yaxis[2] = cosl(beta) * ytemp[2] + sinl(beta) * ztemp[2];
  norm(yaxis);

  Matrix Gamma;
  Gamma.rows = 3;
  Gamma.columns = 3;
  memcpy(Gamma.Element[0], xaxis, 3 * sizeof(double));
  memcpy(Gamma.Element[1], yaxis, 3 * sizeof(double));
  memcpy(Gamma.Element[2], zaxis, 3 * sizeof(double));
  int mode = 0;
  int j = 0;
  int k = 0;
  for (int i = 0; i < 4; i++) {
    mode = 3 * i;
    for (j = 0; j < 3; j++) {
      /*     memcpy(&etran.Element[mode + j][mode], Gamma.Element[j],
                  3 * sizeof(double)); */

      for (k = 0; k < 3; k++) {
        etran.Element[mode + j][mode + k] = Gamma.Element[j][k];
      };
    };
  };

  return etran;
};

SmallMatrix SmallGammaMat(double beta, double xaxis[3]) {
  SmallMatrix etran = {"etran", {0.0}, 12, 12};
  double ztemp[3] = {0.0};
  double ytemp[3] = {0.0};
  double zaxis[3] = {0.0};
  double yaxis[3] = {0.0};
  norm(xaxis);
  int isGimbalLock = (1.0 - fabs(xaxis[1])) <= 1E+06 * tol;
  ytemp[0] =
      -xaxis[1] * isGimbalLock + (-xaxis[0] * xaxis[1]) * (!isGimbalLock);
  ytemp[1] = xaxis[0] * isGimbalLock +
             (xaxis[0] * xaxis[0] + xaxis[2] * xaxis[2]) * (!isGimbalLock);
  ytemp[2] = 0 + (-xaxis[2] * xaxis[1]) * (!isGimbalLock);

  norm(ytemp);

  ztemp[0] =
      (-xaxis[2] * xaxis[0]) * isGimbalLock + -xaxis[2] * (!isGimbalLock);
  ztemp[1] = 0 + (-xaxis[2] * xaxis[1]) * isGimbalLock;
  ztemp[2] = (xaxis[0] * xaxis[0] + xaxis[1] * xaxis[1]) * isGimbalLock +
             xaxis[0] * (!isGimbalLock);
  norm(ztemp);

  zaxis[0] = cosl(beta) * ztemp[0] - sinl(beta) * ytemp[0];
  zaxis[1] = cosl(beta) * ztemp[1] - sinl(beta) * ytemp[1];
  zaxis[2] = cosl(beta) * ztemp[2] - sinl(beta) * ytemp[2];
  norm(zaxis);

  yaxis[0] = cosl(beta) * ytemp[0] + sinl(beta) * ztemp[0];
  yaxis[1] = cosl(beta) * ytemp[1] + sinl(beta) * ztemp[1];
  yaxis[2] = cosl(beta) * ytemp[2] + sinl(beta) * ztemp[2];
  norm(yaxis);

  SmallMatrix Gamma;
  Gamma.rows = 3;
  Gamma.columns = 3;
  memcpy(Gamma.Element[0], xaxis, 3 * sizeof(double));
  memcpy(Gamma.Element[1], yaxis, 3 * sizeof(double));
  memcpy(Gamma.Element[2], zaxis, 3 * sizeof(double));
  int mode = 0;
  int j = 0;
  int k = 0;
  for (int i = 0; i < 4; i++) {
    mode = 3 * i;
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        etran.Element[mode + j][mode + k] = Gamma.Element[j][k];
      };
    };
  };

  return etran;
};

void AssignVariables(StringMatrix *u, int nnodes) {
  char tbf[4];
  for (int i = 0; i < nnodes; i++) {
    sprintf(tbf, "%d", i);
    strcpy(u->Variables[6 * i + 0][0].c, "ux");
    strcat(u->Variables[6 * i + 0][0].c, tbf);
    strcpy(u->Variables[6 * i + 1][0].c, "uy");
    strcat(u->Variables[6 * i + 1][0].c, tbf);
    strcpy(u->Variables[6 * i + 2][0].c, "uz");
    strcat(u->Variables[6 * i + 2][0].c, tbf);

    strcpy(u->Variables[6 * i + 3][0].c, "\u03B8x");
    strcat(u->Variables[6 * i + 3][0].c, tbf);
    strcpy(u->Variables[6 * i + 4][0].c, "\u03B8y");
    strcat(u->Variables[6 * i + 4][0].c, tbf);
    strcpy(u->Variables[6 * i + 5][0].c, "\u03B8z");
    strcat(u->Variables[6 * i + 5][0].c, tbf);
  };
};

void AssignForceVars(StringMatrix *u, int nnodes) {

  char tbf[4];
  for (int i = 0; i < nnodes; i++) {
    sprintf(tbf, "%d", i);
    strcpy(u->Variables[6 * i + 0][0].c, "Fx");
    strcat(u->Variables[6 * i + 0][0].c, tbf);
    strcpy(u->Variables[6 * i + 1][0].c, "Fy");
    strcat(u->Variables[6 * i + 1][0].c, tbf);
    strcpy(u->Variables[6 * i + 2][0].c, "Fz");
    strcat(u->Variables[6 * i + 2][0].c, tbf);

    strcpy(u->Variables[6 * i + 3][0].c, "Mx");
    strcat(u->Variables[6 * i + 3][0].c, tbf);
    strcpy(u->Variables[6 * i + 4][0].c, "My");
    strcat(u->Variables[6 * i + 4][0].c, tbf);
    strcpy(u->Variables[6 * i + 5][0].c, "Mz");
    strcat(u->Variables[6 * i + 5][0].c, tbf);
  };
};

void AssignFixity(const Matrix *K, const Matrix *fixity, const Matrix *concen,
                  Matrix *Kff, StringMatrix *u, Matrix *ends,
                  SmallMatrix RotTrans[MAX_NUMBER_OF_ELEMENTS],
                  SmallMatrix KStiff[MAX_NUMBER_OF_ELEMENTS],
                  int memberID[12][MAX_NUMBER_OF_ELEMENTS]) {

  int ndofs = fixity->rows * fixity->columns;
  int freeIndex[MAX_SIZE]; // map global dof -> free row index
  int isSupport[MAX_SIZE] = {0};
  int nFree = 0;

  Matrix load = {"Fn", {0}, ndofs, 1};
  Matrix Usolution = {"Ufinal", {0.0}, ndofs, 1};

  for (int i = 0; i < fixity->rows; i++) {
    for (int j = 0; j < fixity->columns; j++) {
      int dof = 6 * i + j;

      double fval = fixity->Element[i][j];

      if (isnan(fval)) {
        freeIndex[dof] = nFree;
        nFree++;
      } else {
        freeIndex[dof] = -1;
        Usolution.Element[dof][0] = fval;
        isSupport[i] = 1;
      }

      load.Element[dof][0] = concen->Element[i][j];
    }
  }
  Kff->rows = nFree;
  Kff->columns = nFree;
  StringMatrix uf;
  uf.rows = nFree;
  uf.columns = 1;
  Matrix loadF = {"Fl", {0}, nFree, 1};
  MatrixType m1 = Numerical;

  for (int ig = 0; ig < ndofs; ++ig) {
    int row = freeIndex[ig];
    if (row < 0)
      continue;
    loadF.Element[row][0] = load.Element[ig][0];
    strcpy(uf.Variables[row][0].c, u->Variables[ig][0].c);
    for (int jg = 0; jg < ndofs; ++jg) {
      if (freeIndex[jg] >= 0) {
        int col = freeIndex[jg];
        Kff->Element[row][col] = K->Element[ig][jg];
      }
      loadF.Element[row][0] -= K->Element[ig][jg] * Usolution.Element[jg][0];
    }
  }
  Matrix C;
  C = MatrixVectorSolve(Kff, &uf, &loadF);

  for (int ig = 0; ig < ndofs; ig++) {
    int row = freeIndex[ig];
    if (row >= 0) {
      Usolution.Element[ig][0] = C.Element[freeIndex[ig]][0];
    };
  };
  Matrix R = {"R", {0.0}, ndofs, 1};
  StringMatrix Rname;
  R = Multiply(K, &Usolution);
  AssignForceVars(&Rname, fixity->rows);
  for (int ig = 0; ig < fixity->rows; ig++) {
    if (isSupport[ig]) {
      printf("\nnode %d is a support\n\n", ig);
      for (int jg = 0; jg < 6; jg++) {
        printf("%s = %.18lf\n", Rname.Variables[6 * ig + jg][0].c,
               R.Element[6 * ig + jg][0]);
      };
    };
  };
  printf("\n\n");

  SmallMatrix Tu = {"ugl", {0.0}, 12, 1};
  int il = 0;
  int ol = 0;
  SmallMatrix Fu = {"Flocal", {0.0}, 12, 1};
  double Faxial[MAX_NUMBER_OF_ELEMENTS];
  for (int i = 0; i < ends->rows; i++) {
    for (ol = 0; ol < 12; ol++) {
      Tu.Element[ol][0] = Usolution.Element[memberID[ol][i]][0];
    };
    Tu = SmallMultiply(&RotTrans[i], &Tu);
    Fu = SmallMultiply(&KStiff[i], &Tu);
    Faxial[i] = Fu.Element[0][0];
    fprintf(stdout, "\nMember %d axial force: %.15lf\n", i, Faxial[i]);
  };
  int sos = 0;
  Matrix fluff = {"ufinal", {0.0}, ndofs, 1};
  Matrix Memreact = {"MemberF", {0.0}, ends->rows, 1};
};

Matrix CondenseDOF(const Matrix *K, int *DDOF, int rows, int columns) {
  Matrix Kff = {"Kred", {0.0}, 24, 24};
  /*int ElimDOF[36] = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
                     12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                     26, 27, 28, 32, 33, 34, 38, 39, 40, 44, 45, 46}; */
  int i = 0;
  int j = 0;
  int i1 = 0;
  int j1 = 0;
  int rowdone = 0;

  for (i = 0; i < 24; i++) {
    for (j = 0; j < 24; j++) {
      Kff.Element[i][j] = K->Element[24 + i][24 + j];
    };
  };
  Matrix Kffr = {"Kffr", {0.0}, 20, 20};

  /* for (i = 0; i < 24; i++) {
    for (j = 0; j < 24; j++) {
      if ((i % 6 != 2) && (j % 6 != 2)) {
        Kffr.Element[i1][j1] = Kff.Element[i][j];
        j1++;
        rowdone = 1;
      };
    };
    j1 = 0;
    i1 += (rowdone == 1);
    rowdone = 0;
  }; */

  for (int i = 0; i < 24; ++i) {
    if (i % 6 == 2)
      continue;
    int j1 = 0;
    for (int j = 0; j < 24; ++j) {
      if (j % 6 == 2)
        continue;
      Kffr.Element[i1][j1] = Kff.Element[i][j];
      ++j1;
    }
    ++i1;
  }

  // 24x24
  // 20x20
  // Keeping 12
  // Eliminate 8

  int DDof[12] = {0, 1, 4, 5, 6, 9, 10, 11, 14, 15, 16, 19};
  int ElimDOF[8] = {2, 3, 7, 8, 12, 13, 17, 18};

  Matrix Kee = {"", {0.0}, Kffr.rows - 12, Kffr.columns - 12};

  Matrix Kec = {"", {0.0}, Kffr.rows - 12, 12};

  Matrix Kce = {"", {0.0}, 12, Kffr.columns - 12};

  Matrix Kcc = {"", {0.0}, 12, 12};

  for (i = 0; i < Kee.rows; i++) {
    for (j = 0; j < Kee.columns; j++) {
      Kee.Element[i][j] = Kffr.Element[ElimDOF[i]][ElimDOF[j]];
    };
  };

  for (i = 0; i < Kec.rows; i++) {
    for (j = 0; j < Kec.columns; j++) {
      Kec.Element[i][j] = Kffr.Element[ElimDOF[i]][DDof[j]];
    };
  };

  for (i = 0; i < Kce.rows; i++) {
    for (j = 0; j < Kce.columns; j++) {
      Kce.Element[i][j] = Kffr.Element[DDof[i]][ElimDOF[j]];
    };
  };

  for (i = 0; i < Kcc.rows; i++) {
    for (j = 0; j < Kcc.columns; j++) {
      Kcc.Element[i][j] = Kffr.Element[DDof[i]][DDof[j]];
    };
  };

  Matrix KeeInv = {"", {0.0}, Kee.rows, Kee.columns};
  KeeInv = inverse(&Kee);
  KeeInv = Multiply(&KeeInv, &Kec);
  KeeInv = Multiply(&Kce, &KeeInv);
  Kcc = Subtract(&Kcc, &KeeInv);

  return Kcc;
};

void PrintFreeDOF(int is2d, Matrix *K, Matrix *fixity) {
  int wasprinted = 0;
  int elementsprinted = 0;
  int freeID[MAX_SIZE] = {0};
  int k = 0;
  if (is2d == 1) {
    for (int i = 0; i < fixity->rows; i++) {
      for (int j = 0; j < fixity->columns; j++) {
        if (isnan(fixity->Element[i][j])) {
          freeID[6 * i + j] = 1;
          k++;
        };
      };
    };

    for (int i = 0; i < K->rows; i++) {
      for (int j = 0; j < K->columns; j++) {
        if (freeID[i] && freeID[j]) {
          printf("%+5.5lf\t", K->Element[i][j]);
          wasprinted = 1;
        };
      };
      if (wasprinted == 1) {
        printf("\n");
      };
      wasprinted = 0;
    };
  };

  // printf("\nElements printed = %d\n", elementsprinted);
};

void printvec(double xaxis[3]) {
  printf("<%.6lf  %.6lf  %.6lf>\n", xaxis[0], xaxis[1], xaxis[2]);
};

Matrix AssembleSystemStiffnessMatrix(Matrix *coord_info, Matrix *fixity,
                                     Matrix *properties, Matrix *ends,
                                     StringMatrix *u, Matrix *concen,
                                     Matrix *Gamma1) {
  int nnodes = coord_info->rows;
  int nele = ends->rows;
  if (nnodes * 6 > MAX_SIZE) {
    printf("Matrices too small. Update MAX_SIZE = %d to be at least %d\n",
           MAX_SIZE, nnodes * 6);
    exit(1);
  };

  Matrix K = {"Ksys", {0.0}, nnodes * 6, nnodes * 6};

  // A = [3, 1, 4; 44, 3, 4; 1,4,;];

  if (properties->rows > ends->rows) {
    printf("Fatal error: %d element properties specificed, however ends->rows "
           "= %d ends known\n",
           properties->rows, ends->rows);
    exit(1);
  };

  if (ends->rows > properties->rows) {
    printf("Fatal error: %d members specificed, however only %d properties "
           "specified",
           ends->rows, properties->rows);
    exit(1);
  };

  // A	Izz	  Iyy	J	E	nu	Beta
  int i = 0;
  double xaxis[3];
  int nodeI;
  int nodeJ;
  double A, Izz, Iyy, J, E, v, Beta;
  SmallMatrix RotTrans[MAX_NUMBER_OF_ELEMENTS];
  SmallMatrix RotTransT[MAX_NUMBER_OF_ELEMENTS];
  SmallMatrix LocStiff[MAX_NUMBER_OF_ELEMENTS];
  SmallMatrix KMember[MAX_NUMBER_OF_ELEMENTS];
  for (int hh = 0; hh < MAX_NUMBER_OF_ELEMENTS; hh++) {
    RotTrans[hh] = (SmallMatrix){"RotTrans", {0.0}, 12, 12};
    RotTransT[hh] = (SmallMatrix){"RotTransT", {0.0}, 12, 12};
    LocStiff[hh] = (SmallMatrix){"LocStiff", {0.0}, 12, 12};
    KMember[hh] = (SmallMatrix){"KMember", {0.0}, 12, 12};
  };
  /* SmallMatrix RotTrans = {"RotTrans", {0.0}, 12, 12};
   SmallMatrix RotTransT = {"RotTransT", {0.0}, 12, 12};
   SmallMatrix LocStiff = {"LocStiff", {0.0}, 12, 12};
   SmallMatrix KMember = {"KMember", {0.0}, 12, 12}; */

  int il = 0;
  int ol = 0;
  int memberID[12][MAX_NUMBER_OF_ELEMENTS] = {0};
  double L = 0;
  double As2 = 6.8676;
  double As3 = 18.934;
  MatrixType m1 = Numerical;
  // PRintMatrixData(coord_info);
  for (i = 0; i < nele; i++) {
    nodeI = ends->Element[i][0];
    nodeJ = ends->Element[i][1];
    if (nodeI == nodeJ) {
      printf("Member can't connect to itself");
      exit(1);
    };
    xaxis[0] =
        coord_info->Element[nodeJ - 1][0] - coord_info->Element[nodeI - 1][0];
    xaxis[1] =
        coord_info->Element[nodeJ - 1][1] - coord_info->Element[nodeI - 1][1];
    xaxis[2] =
        coord_info->Element[nodeJ - 1][2] - coord_info->Element[nodeI - 1][2];
    A = properties->Element[i][0];
    Izz = properties->Element[i][1];
    Iyy = properties->Element[i][2];
    J = properties->Element[i][3];
    E = properties->Element[i][4];
    v = properties->Element[i][5];
    Beta = 3.141592653589 * properties->Element[i][6] / 180.0;
    L = sqrt(xaxis[0] * xaxis[0] + xaxis[1] * xaxis[1] + xaxis[2] * xaxis[2]);
    RotTrans[i] = SmallGammaMat(Beta, xaxis);
    // LocStiff[i] = ShearAdjusteElk(A, Izz, Iyy, J, E, v, L, As2, As3);
    LocStiff[i] = Smallelk(A, Izz, Iyy, J, E, v, L);
    RotTransT[i] = SmallTranspose(&RotTrans[i]);
    KMember[i] = SmallMultiply(&LocStiff[i], &RotTrans[i]);
    KMember[i] = SmallMultiply(&RotTransT[i], &KMember[i]);
    for (il = 0; il < 6; il++) {
      memberID[il][i] = (nodeI - 1) * 6 + il;
    };
    for (ol = 6; ol < 12; ol++) {
      memberID[ol][i] = (nodeJ - 1) * 6 + (ol - 6);
    };
    for (ol = 0; ol < 12; ol++) {
      for (il = 0; il < 12; il++) {
        K.Element[memberID[il][i]][memberID[ol][i]] +=
            KMember[i].Element[il][ol];
      };
    };
  };
  Matrix Kff = {"Kff", {0.0}, 0, 0};
  /* Matrix Ks = CondenseDOF(&K, NULL, 12, 12);
   fMatrixPrint(&Ks, NULL, &(MatrixType){Numerical}, stdout);
   Matrix Kgg = Multiply(&Ks, Gamma1);
   Matrix GT = Transpose(Gamma1);
   Kgg = Multiply(&GT, &Kgg);
   fMatrixPrint(&Kgg, NULL, &(MatrixType){Numerical}, stdout); */

  AssignFixity(&K, fixity, concen, &Kff, u, ends, RotTrans, LocStiff, memberID);
  printf("nele = %d\n", nele);

  return K;
};
