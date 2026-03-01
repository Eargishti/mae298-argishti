#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <iostream>
#include <fstream>

#define label size_t
#define scalar long double




struct Point {
scalar x;
scalar y;
scalar z;
};


using Vec3 = std::array<long double, 3>;


template<typename T>
class MeshData {
private:
	std::vector<std::ofstream> Files;

	std::string name;
	std::vector<std::string> st;
	std::vector<std::string> modes;	
public:
	int entries_per_row;
    T *entries;
	size_t N;

	MeshData(T *arg, std::string Names){ 
     this->entries = arg;
	 this->name = Names;
	};
void Open(){
	st = {name + "_binary", name + "_diagnosis", name};
	modes = {"wb", "w", "w"};
	Files.reserve(3);
 	Files[0].open(st.at(0), std::ios::binary); 
for (size_t i = 1; i < st.size(); ++i) {
    const std::string& filename = st[i];
    const char* mode = modes[i].c_str();
	Files[i].open(st[i]);
}

}
void SetEntryAndSizeData(int epr, size_t Number){
	this->entries_per_row = epr;
	this->N = Number;
}

void Write(){
	for (size_t i = 0; i < N; i++){
	 


	}

}

};


class Points{
public:
	FILE *points_binary;
	FILE *points_diagnosis;
	FILE *points;
	std::vector<long double> x;
	std::vector<long double> y;
	std::vector<long double> z;
	int N;
	
void Open(){
	points_binary = fopen("points_binary", "wb");
	points_diagnosis = fopen("points_diagnosis", "w");
	points = fopen("points", "w");
};
void Write(){
	for (int i = 0; i < N; i++){
		fwrite(&x[i], sizeof(long double), 1, points_binary);
		fwrite(&y[i], sizeof(long double), 1, points_binary);
		fwrite(&z[i], sizeof(long double), 1, points_binary);
		fprintf(points_diagnosis, "%.15Lf\t%.15Lf\n", x[i], y[i]);
		fprintf(points, "( %.15Lf %.15Lf %.15Lf )\n", x[i], y[i], z[i]);
	};
};
};

class Faces{
public:
	FILE *faces_binary;
	FILE *faces_diagnosis;
	FILE *faces;
	int N;
	int **p;
	int vertices = 6;
void Open(){
	faces_binary = fopen("faces_binary", "wb");
	faces_diagnosis = fopen("faces_diagnosis", "w");
	faces = fopen("faces", "w");
};
private:
void Write(){
	for (int i = 0; i < N; i++){
		fwrite(&vertices, sizeof(int), 1, faces_binary);
		fprintf(faces_diagnosis, "%d\t", vertices);
		fprintf(faces, "%d(", vertices);
         for (int ii = 0; ii < vertices; ii++){
			fwrite(&p[i][ii], sizeof(int), 1, faces_binary);
			fwrite(&p[i][ii], sizeof(int), 1, faces_binary);
			fprintf(faces_diagnosis, "%d\t", p[i][ii]);
			fprintf(faces, "%d ", p[i][ii]);
		
		 };
		fprintf(faces_diagnosis, "\n");
		fprintf(faces, ")\n");
	};
};

public:
void TriangleWrite(int num){
	N = num;
	vertices = 3;
	Write();
};
void QuadrilateralWrite(int num){
	N = num;
	vertices = 4;
	Write();
};
};


class Cells {
private:
	FILE *cells_binary;
	FILE *cells_diagnosis;
	FILE *cells;
	int faces_per_cell = 6;
	int **faces;
	
public:
	int N;
void Open(){
	cells_binary = fopen("cells_binary", "wb");
	cells_diagnosis = fopen("cells_diagnosis", "w");
	cells = fopen("cells", "w");
};
void Write(){
	for (int i = 0; i < N; i++){
		fwrite(&faces_per_cell, sizeof(int), 1, cells_binary);
		fprintf(cells_diagnosis, "%d\t", faces_per_cell);
		fprintf(cells, "%d(", faces_per_cell);
         for (int ii =0; ii < faces_per_cell; ii++){
			fwrite(&faces[i][ii], sizeof(int), 1, cells_binary);
			fprintf(cells_diagnosis, "%d\t", faces[i][ii]);
			fprintf(cells, "%d ", faces[i][ii]);
		
		 };
		fprintf(cells_diagnosis, "\n");
		fprintf(cells, ")\n");
	};

};
void HexahedronWrite(){
	faces_per_cell = 6;
	Write();
};
void WedgeWrite(){
	faces_per_cell = 5;
	Write();
};
};


class Owner{
	FILE *owner_binary;
	FILE *owner_diagnosis;
	FILE *owner;
	int *clabel;
	int N;
	
void Open(){
	owner_binary = fopen("owner_binary", "wb");
	owner_diagnosis = fopen("owner_diagnosis", "w");
	owner = fopen("owner", "w");
}
void Write(int *clabel, int N){
		fwrite(clabel, sizeof(int), N, owner_binary);
	for (int i = 0; i < N; i++){
		fprintf(owner_diagnosis, "%d\n", clabel[i]);
		fprintf(owner, "%d\n", clabel[i]);
	}

}
};


class Neighbour{
private:
	FILE *neighbour_binary;
	FILE *neighbour_diagnosis;
	FILE *neighbour;
	int *clabel;
	int N;
public:
void Open(){
	neighbour_binary = fopen("neighbour_binary", "wb");
	neighbour_diagnosis = fopen("neighbour_diagnosis", "w");
	neighbour = fopen("neighbour", "w");
}
void Write(int *clabel, int N){
		fwrite(clabel, sizeof(int), N, neighbour_binary);
	for (int i = 0; i < N; i++){
		fprintf(neighbour_diagnosis, "%d\n", clabel[i]);
		fprintf(neighbour, "%d\n", clabel[i]);
	}
}
};
