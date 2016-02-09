#include <iostream>
#include <iomanip>

using namespace std;

void displayField(int row, int col, float* field);
void displayVectorField(int row, int col, float* u, float* v);
void computeCurls_uniform(float * w, float * u, float * v);