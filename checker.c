#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef enum { North1, North2, West1, West2, IN1, OUT1 } fcode;

typedef enum { AEROFOIL, NORTH, SOUTH, EAST, WEST, IN, OUT, INLET, OUTLET, TOP, BOTTOM, REVERSE, TOP2, BOTTOM2 } facecode;

typedef enum { FrontAndBack, Aerofoil, Inlet, Outlet, Top1, Bottom1, Top2, Bottom2 } BoundaryPatch;

typedef enum { value1, value2 } ValueType;

typedef struct CodeMap {

  ValueType val1;
  union {
    facecode f1;
    BoundaryPatch b1;

  } Map1;

} CodeMap;

const char *stringfacecode(CodeMap c1) {
  switch (c1.val1) {
  case 0:
    switch (c1.Map1.f1) {
    case AEROFOIL:
      return "AEROFOIL";
    case NORTH:
      return "NORTH";
    case SOUTH:
      return "SOUTH";
    case EAST:
      return "EAST";
    case WEST:
      return "WEST";
    case IN:
      return "IN";
    case OUT:
      return "OUT";
    case INLET:
      return "INLET";
    case OUTLET:
      return "OUTLET";
    case TOP:
      return "TOP";
    case BOTTOM:
      return "BOTTOM";
    case REVERSE:
      return "REVERSE";
    case TOP2:
      return "TOP2";
    case BOTTOM2:
      return "BOTTOM2";
    };
    break;
  case 1:
    switch (c1.Map1.b1) {
    case FrontAndBack:
      return "FrontAndBack";
    case Aerofoil:
      return "Aerofoil";
    case Inlet:
      return "Inlet";
    case Outlet:
      return "Outlet";
    case Top1:
      return "Top1";
    case Bottom1:
      return "Bottom1";
    case Top2:
      return "Top2";
    case Bottom2:
      return "Bottom2";
    };
    break;
  };
};

void getcell(int ID, FILE *cells, int *f) {
  if (cells == NULL) {
    printf("failed to open file\n");
  }

  long long off = (long long)(ID) * 7 * sizeof(int);

  if (fseek(cells, off, SEEK_SET) != 0) {
    printf("fseek failed to run\n");
    return;
  }

  if (fread(f, sizeof(int), 7, cells) != 7) {
    printf("\nfread failed to run\n");
    return;
  }
}

typedef struct Point {
  long double x;
  long double y;
  long double z;
  int plabel;

} Point;

typedef struct Face {
  Point points[4];
  facecode bcode;
  int owner;
  int neighbour;

} Face;

typedef struct Cell {
  Face faces[6];
  int clabel;

} Cell;

void facefetch(int *flabel, FILE *faces, FILE *neighbour, FILE *owner, FILE *cells, int *f, int *nFaces, int *startFace, int p[4]) {
  // Per-face binary layout in "faces": 4 vertex ints + 1 bcode + 1 faceID (=6 ints)
  const long long FACE_REC_SIZE = 6LL * sizeof(int);
  // Per-face binary layout in "owner": 1 clabel + 1 bcode (=2 ints)
  const long long OWNER_REC_SIZE = 2LL * sizeof(int);

  if (!faces || !owner || !cells) {
    printf("null file pointer(s)\n");
    return;
  }

  int faceID = 0;
  int bcode = 0;
  int clabel = 0; // owner cell id from owner file
  CodeMap map1;
  map1.val1 = 0;

  // --- seek/read faces record ---
  long long off = (long long)(*flabel) * FACE_REC_SIZE;
  if (fseek(faces, off, SEEK_SET) != 0) {
    perror("fseek faces failed");
    return;
  }
  if (fread(p, sizeof(int), 4, faces) != 4) {
    printf("Face label %d failed to read 4 points\n", *flabel);
    return;
  }
  if (fread(&bcode, sizeof(int), 1, faces) != 1) {
    printf("Face label %d failed to read bcode\n", *flabel);
    return;
  }
  if (fread(&faceID, sizeof(int), 1, faces) != 1) {
    printf("Face label %d failed to read faceID\n", *flabel);
    return;
  }

  // --- seek/read owner record ---
  off = (long long)(*flabel) * OWNER_REC_SIZE;
  if (fseek(owner, off, SEEK_SET) != 0) {
    perror("fseek owner failed");
    return;
  }
  if (fread(&clabel, sizeof(int), 1, owner) != 1) {
    printf("Face label %d failed to read owner cell label\n", *flabel);
    return;
  }
  int bc_owner = 0;
  if (fread(&bc_owner, sizeof(int), 1, owner) != 1) {
    printf("Face label %d failed to read owner bcode\n", *flabel);
    return;
  }

  getcell(clabel, cells, f);
  map1.Map1.f1 = bcode;
  printf("Face %d has an owner of %d (faceID=%d, bcode=%d, face=%s)\n", *flabel, clabel, faceID, bcode, stringfacecode(map1));
  printf("Face %d: %d %d %d %d \n", *flabel, p[0], p[1], p[2], p[3]);
  printf("Cell %d:\n", f[6]);
  printf("%d\t%d\t%d\t%d\t%d\t%d\n", f[0], f[1], f[2], f[3], f[4], f[5]);
  if (*flabel < nFaces[0] + startFace[0]) {
    // FrontAndBack example logic…
  }
}

void getboundarydomain(FILE *boundary, int nFaces[], int startFace[], int JMax) {
  char stong[4][32];
  for (int J = 0; J < JMax; J++) {

    fscanf(boundary, "%s\t%d\t%d\n", stong[J], &nFaces[J], &startFace[J]);
    printf("%s\t%d\t%d\n", stong[J], nFaces[J], startFace[J]);
  };
}

int main() {
  FILE *boundary;
  FILE *cells;
  FILE *faces;
  FILE *owner;
  FILE *neighbour;
  cells = fopen("cellcheck", "rb");
  boundary = fopen("Truebound", "r");
  faces = fopen("facecheck", "rb");
  owner = fopen("owncheck", "rb");
  if (cells == NULL || faces == NULL || owner == NULL || boundary == NULL) {
    printf("\nFile failed to open\n");
    return 3;
  };

  if (!boundary) {
    printf("File failed to open\n");
    return 4;
  };
  int JMax = 4;
  int *nFaces = (int *)malloc(JMax * sizeof(int));
  int *startFace = (int *)malloc(JMax * sizeof(int));
  int f[7];
  int flabel;
  getboundarydomain(boundary, nFaces, startFace, JMax);
  printf("\nSelect a face to check: ");
  scanf("%d", &flabel);
  int p[4];
  facefetch(&flabel, faces, neighbour, owner, cells, f, nFaces, startFace, p);
};
