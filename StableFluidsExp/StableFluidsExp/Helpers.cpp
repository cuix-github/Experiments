#include "Helpers.h"

void displayField(int row, int col, float* field){
	cout << std::setprecision(4) << setiosflags(ios::fixed);
	for (int i = 0; i != row; i++){
		for (int j = 0; j != col; j++){
			cout << std::setw(7) << field[i * col + j] << ", ";
			if (j == col - 1)
				cout << endl;
		}
	}
}

void displayVectorField(int row, int col, float* u, float* v){
	cout << std::setprecision(4) << setiosflags(ios::fixed);
	for (int i = 0; i != row; i++){
		for (int j = 0; j != col; j++){
			cout << std::setw(3) << "v(" << u[i * col + j] << "," << v[i * col + j] << ") ";
			if (j == col - 1)
				cout << endl;
		}
	}
}

void displayFieldInv(int row, int col, float * field){
	cout << std::setprecision(4) << setiosflags(ios::fixed);
	for (int i = row - 1; i >= 0; i--){
		for (int j = 0; j != col; j++){
			cout << std::setw(7) << field[i * col + j] << ", ";
			if (j == col - 1)
				cout << endl;
		}
	}
}

void displayVectorFieldInv(int row, int col, float * u, float * v){
	cout << std::setprecision(4) << setiosflags(ios::fixed);
	for (int i = row - 1; i >= 0; i--){
		for (int j = 0; j != col; j++){
			cout << std::setw(3) << "v(" << u[i * col + j] << "," << v[i * col + j] << ") ";
			if (j == col - 1)
				cout << endl;
		}
	}
}

void computeCurls_uniform(int dim, float * w, float * u, float * v)
{
	float h = 1.0f / dim;

	int N = dim - 2;

	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			float du, dv, curl;
			float uiplus1, uiminus1;
			float u0, v0;
			u0 = u[IX(i, j)];
			v0 = v[IX(i, j)];
			uiplus1 = u[IX(i + 1, j)];
			uiminus1 = u[IX(i - 1, j)];
			du = 0.5f * (u[IX(i + 1, j)] - u[IX(i - 1, j)]) / h;
			dv = 0.5f * (v[IX(i, j + 1)] - v[IX(i, j - 1)]) / h;
			curl = dv - du;
			w[IX(i, j)] = curl;
		}
	}
}

void computeDivergence_unifrom(int N, float * u, float * v, float * div){

	for (int i = 1; i <= N; i++){
		for (int j = 1; j <= N; j++){
			float du, dv;
			float divergence;
			du = -0.5f * (u[IX(i + 1, j)] - u[IX(i - 1, j)]) / N;
			dv = -0.5f * (v[IX(i, j + 1)] - v[IX(i, j - 1)]) / N;
			divergence = du + dv;
			div[IX(i, j)] = divergence;
		}
	}
}