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

void displayField(int row, int col, float * field, float * fieldOrigin)
{
	cout << endl << "Current timestep field:" << endl;
	displayField(row, col, field);
	cout << endl << "Previous timestep field:" << endl;
	displayField(row, col, fieldOrigin);
}

void displayVectorField(int row, int col, float* u, float* v){
	cout << std::setprecision(4) << setiosflags(ios::fixed);
	int N = row;
	for (int i = 0; i != row; i++){
		for (int j = 0; j != col; j++){
			cout << std::setw(3) << "v(" << u[i * col + j] << "," << v[i * col + j] << ") ";
			if (j == col - 1)
				cout << endl;
		}
	}
}

void displayVectorField(int row, int col, float * u, float * v, float * u0, float * v0)
{
	cout << endl << "Current timestep field:" << endl;
	displayVectorField(row, col, u, v);
	cout << endl << "Previous timepstep field:" << endl;
	displayVectorField(row, col, u0, v0);
}

void computeCurls_uniform(int N, float * w, float * u, float * v)
{
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			float du, dv;
			float uiplus1, uiminus1;
			float u0, v0;
			du = 0.5f * (u[IX(i, j + 1)] - u[IX(i, j - 1)]) * N;
			dv = 0.5f * (v[IX(i + 1, j)] - v[IX(i - 1, j)]) * N;
			w[IX(i, j)] = dv - du;
		}
	}
}

void vector_potential_inv_2D(int N, float * u, float * v, float * psi)
{
	float dpsi_dy, dpsi_dx;
	
	// wx and wy finally will be used for updating velocity field (u = wx, v = wy)
	// While 3D Vector field is constructed as (u = 0, v = 0, w = psi)
	// Where psi is the stream function computed from the linear solver.

	// TODO: 2 Nested loop for wx
	for (int i = 1; i <= N; i++){
		for (int j = 1; j <= N; j++){
			dpsi_dy = 0.5f * (psi[IX(i + 1, j)] - psi[IX(i - 1, j)]) * N;
			u[IX(i, j)] = -dpsi_dy;
		}
	}

	// TODO: 2 Nested loop for wy
	for (int i = 1; i <= N; i++){
		for (int j = 1; j <= N; j++){
			dpsi_dx = 0.5f * (psi[IX(i, j + 1)] - psi[IX(i, j - 1)]) * N;
			v[IX(i, j)] = dpsi_dx;
		}
	}
}

void linear_combine_sub(int N, float * f_out, float * f, float * f0)
{
	for (int i = 1; i <= N; i++){
		for (int j = 1; j <= N; j++){
			f_out[IX(i, j)] = f[IX(i, j)] - f0[IX(i, j)];
		}
	}
}

void linear_combine_add(int N, float * f_out, float * f, float * f0)
{
	for (int i = 1; i <= N; i++){
		for (int j = 1; j <= N; j++){
			f_out[IX(i, j)] = f[IX(i, j)] - f0[IX(i, j)];
		}
	}
}

void zeros(int N, float * field){
	for (int i = 0; i != N + 2; i++){
		for (int j = 0; j != N + 2; j++){
			field[IX(i, j)] = 0.0f;
		}
	}
}

void scaler(int N, float *field, float factor)
{
	for (int i = 1; i <= N; i++){
		for (int j = 1; j <= N; j++){
			field[IX(i, j)] *= factor;
		}
	}
}

double relative_error(int N, float * curr, float * prev)
{
	double error = 0.0;
	int n = 0;

	for (int i = 1; i <= N; i++){
		for (int j = 1; j <= N; j++){
			if (curr[IX(i, j)] != 0){
				if (curr[IX(i, j)] != prev[IX(i, j)])
				{
					error += std::abs(1 - prev[IX(i, j)] / curr[IX(i, j)]);
					++n;
				}
			}
		}
	}

	if (n != 0)
		error /= n;
	return error;
}

void computeDivergence_unifrom(int N, float * u, float * v, float * div){
	int i, j;
	FOR_EACH_CELL
		div[IX(i, j)] = - 0.5f*(u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]) / N;
	END_FOR
}