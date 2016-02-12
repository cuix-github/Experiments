#include <iostream>
#include <iomanip>

using namespace std;

#define IX(i,j) ((i)+(N+2)*(j))

void displayField(int row, int col, float* field);
void displayVectorField(int row, int col, float* u, float* v);
void computeCurls_uniform(int N, float h, float * w, float * u, float * v);