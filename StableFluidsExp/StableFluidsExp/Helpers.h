#include <iostream>
#include <iomanip>
#include <math.h>

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

void displayField(int row, int col, float * field);
void displayField(int row, int col, float * field, float * fieldOrigin);
void displayVectorField(int row, int col, float * u, float * v);
void displayVectorField(int row, int col, float * u, float * v, float * u0, float * v0);
void zeros(int N, float * field);
void scaler(int N, float *field, float factor);
double relative_error(int N, float * curr, float * prev);
void computeCurls_uniform(int N, float * w, float * u, float * v);
void find_vector_potential_2D(int N, float * u, float * v, float * psi);
void linear_combine_add(int N, float * f_out, float * f, float * f0);
void linear_combine_sub(int N, float * f_out, float * f, float * f0);
void computeDivergence_unifrom(int N, float * u, float * v, float * div);
float lerp(float t, float x0, float x1);
float bilerp(float t, float s, float x0, float x1, float y0, float y1);

typedef struct {
	float x, y;
	float padding;
}Particle;