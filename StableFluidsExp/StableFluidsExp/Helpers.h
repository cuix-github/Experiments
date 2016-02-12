#include <iostream>
#include <iomanip>

using namespace std;

#define IXCURL(i,j) ((i) * (N + 3) + (j))
#define IX(i,j) ((j)+(N+2)*(i))

void displayField(int row, int col, float* field);
void displayVectorField(int row, int col, float* u, float* v);
void displayVectorFieldInv(int row, int col, float * u, float * v);
void displayFieldInv(int row, int col, float * field);
void computeCurls_uniform(int N, float * w, float * u, float * v);