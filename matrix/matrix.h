#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define MAX_SIZE 50
#define MAX_NUMBER_OF_MATRICES 40
#define BUFFER_SIZE 256
#define tol 1E-13
#define VARIABLE_SIZE 8
#define MAX_NUMBER_OF_ELEMENTS 10

typedef struct {

  char name[BUFFER_SIZE];
  double Element[MAX_SIZE][MAX_SIZE];
  uint8_t rows;
  uint8_t columns;

} Matrix;

typedef struct {
char name[BUFFER_SIZE];
double Element[12][12];
uint8_t rows;
uint8_t columns;


} SmallMatrix;


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

SmallMatrix SmallMultiply(const SmallMatrix *A, const SmallMatrix *B);

SmallMatrix SmallTranspose(const SmallMatrix *matrix);

void SymmetricSmallUT(SmallMatrix *matrix);

void infoprint();

void AllocateFiles(int argc, char *argv[], FILE ***matrixfiles);


void fMatrixPrint(const Matrix *matrices, StringMatrix *stringmatrix, MatrixType *m1,FILE *matrixfile);;

void PrintMatrixData(Matrix *matrices, StringMatrix *stringmatrix, MatrixType *m1);


void SidebySide(Matrix *matrix1, Matrix *matrix2, FILE *f1);
void PRintMatrixData(const Matrix *matrices);


void PrintFileContent(FILE **matrixfiles, int filecount);

int IsLetter(char a);

int IsNumber(char a);

int IsWhiteSpace(char a);

void GetRowData(char *Buffer, Matrix *matrix1, FILE **matrixfile, int row, mode1 mode, MatrixType *m1, StringMatrix *stringmat);

void GetVariableName(char *string1, Matrix *matrix1, char *filename, int *linenumber, FILE *matrixfile, MatrixType *m1, StringMatrix *stringmat1);

void GetToType(FILE *matrixfile, Matrix *matrices, int MatrixID, char *Buffer, StringMatrix *stringmat);

void SaveFileMatrixData(FILE *matrixfile, Matrix *matrices, int *MatrixID, char *filename, StringMatrix *stringmats, int *count, int m[MAX_SIZE]);

Matrix Multiply(const Matrix *A,const Matrix *B);

double Det(const Matrix *matrix, Matrix *U, int printflag, Matrix *colvector, int Swapper[MAX_SIZE][2], int *swapcount);

double SimpleDet(const Matrix *matrix, Matrix *U, int printflag);

double Tran(Matrix matrix, Matrix *T);

Matrix Transpose(const Matrix *matrix);

double Inverse(const Matrix *matrix, Matrix *Inv);

Matrix inverse(const Matrix *matrix);

void Add(Matrix *A, Matrix *B);

Matrix Subtract(const Matrix *A,const Matrix *B);

void CreateRandomMatrix(Matrix *matrix, unsigned int rows, unsigned int columns, int max, int scale, int *ID);

void InitializeValues(Matrix *matrix);

Matrix MatrixVectorSolve(Matrix *matrix, StringMatrix *vars, Matrix *cols);

void MatrixEnumForm(FILE *header, Matrix *matrices, StringMatrix *stringmats, int *ID, int m[MAX_SIZE]);

//Homework
Matrix elk(double A, double Izz, double Iyy, double J, double E, double nu, double L);

Matrix GammaMat(double beta, double xaxis[3]);


SmallMatrix Smallelk(double A, double Izz, double Iyy, double J, double E, double nu, double L);

SmallMatrix SmallGammaMat(double beta, double xaxis[3]);


//Project rows by columns,
//concen, fixity, 
void AssignFixity(const Matrix *K, const Matrix *fixity, const Matrix *concen, Matrix *Kff, StringMatrix *u, Matrix *ends, SmallMatrix RotTrans[MAX_NUMBER_OF_ELEMENTS], SmallMatrix KStiff[MAX_NUMBER_OF_ELEMENTS], int memberID[12][MAX_NUMBER_OF_ELEMENTS]);

void PrintFreeDOF(int is2D, Matrix *K, Matrix *fixity);

void MatrixDefineForm(FILE *header, Matrix *matrices, StringMatrix *stringmats, int *ID, int m[MAX_SIZE]);

void AssignVariables(StringMatrix *u, int nnodes);

Matrix AssembleSystemStiffnessMatrix(Matrix *coord_info, Matrix *fixity, Matrix *properties, Matrix *ends, StringMatrix *u, Matrix *concen);

void AssignForceVars(StringMatrix *u, int nnodes);

Matrix DebugMultiply(const Matrix *A, const Matrix *B, int problematicrow[2]);
