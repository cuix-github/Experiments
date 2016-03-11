#include <iostream>
#include <iomanip>

using namespace std;

#define IXCURL(i,j) ((i) * (N + 2) + (j))
#define IX(i,j) ((i)*(N+2)+(j))
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
#define END_FOR }}

#ifdef __cplusplus__
#include <cstdlib>
#else
#include <stdlib.h>
#endif

void displayField(int row, int col, float* field);
void displayVectorField(int row, int col, float* u, float* v);
void displayVectorFieldInv(int row, int col, float * u, float * v);
void displayFieldInv(int row, int col, float * field);
void zeros(int N, float * field);
void computeCurls_uniform(int N, float * w, float * u, float * v);
void computeDivergence_unifrom(int N, float * u, float * v, float * div);