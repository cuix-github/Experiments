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

float lerp(float t, float x0, float x1){
	return (1 - t) * x0 + t * x1;
}

float bilerp(float s, float t, float x0, float x1, float y0, float y1){
	return lerp(t, lerp(s, x0, x1), lerp(s, y0, y1));
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
			du = 0.5f * (u[IX(i + 1, j)] - u[IX(i - 1, j)]) * N;
			dv = 0.5f * (v[IX(i, j + 1)] - v[IX(i, j - 1)]) * N;
			w[IX(i, j)] = dv - du;
		}
	}
}

void computeCurls_uniform_smoother(int N, float * w, float * u, float * v)
{
	for (int i = 1; i <= N; i++){
		for (int j = 1; j <= N; j++){
			float average_x_top = 0.5 * (u[IX(i, j)] + u[IX(i, j + 1)]);
			float average_x_bottom = 0.5 * (u[IX(i + 1, j)] + u[IX(i, j + 1)]);
			float average_y_left = 0.5 * (v[IX(i, j)] + v[IX(i + 1, j)]);
			float average_y_right = 0.5 * (v[IX(i, j + 1)] + v[IX(i + 1, j + 1)]);
			float du = average_x_bottom - average_x_top;
			float dv = average_y_right - average_y_left;
			w[IX(i, j)] = dv - du;
		}
	}
}

void find_vector_potential_2D(int N, float * u, float * v, float * psi)
{
	float dpsi_dy, dpsi_dx;
	zeros(N, u);
	zeros(N, v);
	
	// wx and wy finally will be used for updating velocity field (u = wx, v = wy)
	// While 3D Vector field is constructed as (u = 0, v = 0, w = psi)
	// Where psi is the stream function computed from the linear solver.

	// TODO: 2 Nested loop for wx
	for (int i = 1; i <= N; i++){
		for (int j = 1; j <= N; j++){
			dpsi_dy = 0.5f * (psi[IX(i + 1, j)] - psi[IX(i - 1, j)]) * N;
			dpsi_dx = 0.5f * (psi[IX(i, j + 1)] - psi[IX(i, j - 1)]) * N;
			u[IX(i, j)] = dpsi_dy;
			v[IX(i, j)] = -dpsi_dx;
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
			f_out[IX(i, j)] = f[IX(i, j)] + f0[IX(i, j)];
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
		return error /= n;
}

void computeDivergence_unifrom(int N, float * u, float * v, float * div){
	int i, j;
	FOR_EACH_CELL
		div[IX(i, j)] = 0.5f * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]) / N;
	END_FOR
}

float interpolate(int N, float x, float y, float * field){
	int i, j;
	float fx, fy;

	get_barycentric(x, i, fx, 0, N);
	get_barycentric(y, j, fy, 0, N);

	return bilerp(fx, fy, field[IX(i, j)], field[IX(i + 1, j)], field[IX(i, j + 1)], field[IX(i + 1, j + 1)]);
}