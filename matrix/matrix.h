#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define MAX_SIZE 50
#define MAX_NUMBER_OF_MATRICES 50
#define BUFFER_SIZE 256
#define tol 1E-13
#define VARIABLE_SIZE 8

typedef struct {

  char name[BUFFER_SIZE];
  double Element[MAX_SIZE][MAX_SIZE];
  uint8_t rows;
  uint8_t columns;

} Matrix;


typedef enum { variable_name, matrix_type, matrix_values } Reading_Type;

typedef enum { Brackets, NoBrackets } mode1;

typedef enum { Numerical, Symbolic } MatrixType;

typedef struct {
  char c[VARIABLE_SIZE + 1];
} char3;

typedef struct {
  char name[BUFFER_SIZE];
  char3 Variables[MAX_SIZE][MAX_SIZE];
  uint8_t rows;
  uint8_t columns;

} StringMatrix;


void infoprint();

void AllocateFiles(int argc, char *argv[], FILE ***matrixfiles);


void fMatrixPrint(Matrix *matrices, StringMatrix *stringmatrix, MatrixType *m1,FILE *matrixfile);;

void PrintMatrixData(Matrix *matrices, StringMatrix *stringmatrix, MatrixType *m1);


void SidebySide(Matrix *matrix1, Matrix *matrix2, FILE *f1);
void PRintMatrixData(Matrix *matrices, int k);


void PrintFileContent(FILE **matrixfiles, int filecount);

int IsLetter(char a);

int IsNumber(char a);

int IsWhiteSpace(char a);

void GetRowData(char *Buffer, Matrix *matrix1, FILE **matrixfile, int row, mode1 mode, MatrixType *m1, StringMatrix *stringmat);

void GetVariableName(char *string1, Matrix *matrix1, char *filename, int *linenumber, FILE *matrixfile, MatrixType *m1, StringMatrix *stringmat1);

void GetToType(FILE *matrixfile, Matrix *matrices, int MatrixID, char *Buffer, StringMatrix *stringmat);

void SaveFileMatrixData(FILE *matrixfile, Matrix *matrices, int *MatrixID, char *filename, StringMatrix *stringmats, int *count, int m[MAX_SIZE]);

void Multiply(Matrix *A, Matrix *B, Matrix *D);

double Det(const Matrix *matrix, Matrix *U, int printflag, Matrix *colvector, int Swapper[MAX_SIZE][2], int *swapcount);

double SimpleDet(const Matrix *matrix, Matrix *U, int printflag);

double Tran(Matrix matrix, Matrix *T);

double Inverse(Matrix *matrix, Matrix *Inv);

void Add(Matrix *A, Matrix *B);

void Subtract(Matrix *A, Matrix *B);

void CreateRandomMatrix(Matrix *matrix, unsigned int rows, unsigned int columns, int max, int scale, int *ID);

void InitializeValues(Matrix *matrix);

void MatrixVectorSolve(Matrix *matrix, StringMatrix *vars, Matrix *cols);

void MatrixEnumForm(FILE *header, Matrix *matrices, StringMatrix *stringmats, int *ID, int m[MAX_SIZE]);

//Homework
Matrix elk(long double A, long double Izz, long double Iyy, long double J, long double E, long double nu, long double L);

Matrix GammaMat(long double beta, long double xaxis[3]);

//Project rows by columns,

void AssignFixity();

