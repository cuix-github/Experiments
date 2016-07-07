#include "helpers.h"

void displayScalarField(int row, int col, float* field){
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

void computeCurls_uniform(int N, float * w, float * u, float * v)
{
	float h = 1.0f / N;
	//for (int i = 1; i <= N; i++)
	//{
	//	for (int j = 1; j <= N; j++)
	//	{	
	//		// Stokes Theorem
	//		float coef = 1.0f / (8 * (std::pow(h, 2)));
	//		float du1 = h * (u[IX(i - 1, j - 1)] + 2 * u[IX(i, j - 1)] + u[IX(i + 1, j - 1)]);
	//		float du0 = h * (u[IX(i + 1, j + 1)] + 2 * u[IX(i, j + 1)] + u[IX(i - 1, j + 1)]);
	//		float dv1 = h * (v[IX(i + 1, j - 1)] + 2 * v[IX(i + 1, j)] + v[IX(i + 1, j + 1)]);
	//		float dv0 = h * (v[IX(i - 1, j + 1)] + 2 * v[IX(i - 1, j)] + v[IX(i - 1, j - 1)]);
	//		w[IX(i, j)] = coef * (du1 + dv1 - du0 - dv0);
	//	}
	//}

	//// Least Square
	//for (int i = 2; i <= N - 1; i++)
	//{
	//	for (int j = 2; j <= N - 1; j++)
	//	{
	//		float coef = 1 / (10 * h);
	//		float du = 2 * u[IX(i, j + 2)] + u[IX(i, j + 1)] - u[IX(i, j - 1)] - 2 * u[IX(i, j - 2)];
	//		float dv = 2 * v[IX(i + 2, j)] + v[IX(i + 1, j)] - v[IX(i - 1, j)] - 2 * v[IX(i - 2, j)];
	//		w[IX(i, j)] = coef * dv - coef * du;
	//	}
	//}

	// Richardson Extrapolation
	for (int i = 2; i <= N - 1; i++)
	{
		for (int j = 2; j <= N - 1; j++)
		{
			float coef = 1 / (12 * h);
			float du = -v[IX(i, j + 2)] + 8 * v[IX(i, j + 1)] - 8 * v[IX(i, j - 1)] - v[IX(i, j - 2)];
			float dv = -u[IX(i + 2, j)] + 8 * u[IX(i + 1, j)] - 8 * u[IX(i - 1, j)] - u[IX(i - 2, j)];
			w[IX(i, j)] = coef * dv - coef * du;
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

void computeDivergence_unifrom(int N, float * u, float * v, float * div){
	LOOP_CELLS{
		div[IX(i, j)] = 0.5f * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]) / N;
	}
}

void computeBuoyancy(int N, float * v, float * d, float * ambientTemperature, float alpha, float beta, float dt)
{
	float sampledDensity, sampledTemperature;
	LOOP_CELLS{
		sampledDensity = d[IX(i, j)];
		sampledTemperature = ambientTemperature[IX(i, j)];
		v[IX(i, j)] = v[IX(i, j)] - alpha * sampledDensity + beta * sampledTemperature;
	}
}
void computeVortConf(int N, float * u, float * v, float dt, float vort_conf_eps)
{
	float * vort = (float*)malloc((N + 2) * (N + 2) * sizeof(float));
	float * fvortu = (float*)malloc((N + 2) * (N + 2) * sizeof(float));
	float * fvortv = (float*)malloc((N + 2) * (N + 2) * sizeof(float));
	if (!vort || !fvortu || !fvortv)
		exit(0);

	zeros(N, vort);
	zeros(N, fvortu);
	zeros(N, fvortv);

	vec3 dir_vort = vec3(0.0f, 0.0f, 0.0f);
	vec3 vorticity = vec3(0.0f, 0.0f, 0.0f);
	
	// Get curls
	computeCurls_uniform(N, vort, u, v);

	LOOP_CELLS {
		float gradu, gradv, gradlen;
		gradu = 0.5f * N * (vort[IX(i + 1, j)] - vort[IX(i - 1, j)]);
		gradv = 0.5f * N * (vort[IX(i, j + 1)] - vort[IX(i, j - 1)]);
		gradlen = std::sqrt(std::pow(gradu, 2) + std::pow(gradv, 2));
		gradu *= 1.0f / (gradlen + 10e-20);
		gradv *= 1.0f / (gradlen + 10e-20);
		dir_vort.x = gradu;
		dir_vort.y = gradv;
		dir_vort.z = 0.0f;
		vorticity.x = vorticity.y = 0.0f;
		vorticity.z = vort[IX(i, j)];
		vec3 fconf = vec3(0.0f, 0.0f, 0.0f);
		fconf = cross(vorticity, dir_vort);
		fconf.x *= (vort_conf_eps * (1.0f / N));
		fconf.y *= (vort_conf_eps * (1.0f / N));
		fconf.z *= (vort_conf_eps * (1.0f / N));
		
		if (fconf.z > 0 || fconf.z < 0)
		{
			cout << endl << "Vorticity confinement calculation errors" << endl;
			system("Pause");
			exit(0);
		}

		fvortu[IX(i, j)] = fconf.x;
		fvortv[IX(i, j)] = fconf.y;
	}

	LOOP_CELLS {
		u[IX(i, j)] += (fvortu[IX(i, j)] + fvortu[IX(i - 1, j)]) * 0.5f;
		v[IX(i, j)] += (fvortv[IX(i, j)] + fvortv[IX(i, j - 1)]) * 0.5f;
	}

	free(vort);
}


float interpolate(int N, float x, float y, float * field){
	int i, j;
	float fx, fy;

	get_barycentric(x, i, fx, 0, N);
	get_barycentric(y, j, fy, 0, N);

	return bilerp(fx, fy, field[IX(i, j)], field[IX(i + 1, j)], field[IX(i, j + 1)], field[IX(i + 1, j + 1)]);
}

vec3 cross(vec3 v1, vec3 v2){
	vec3 result = vec3(0.0f, 0.0f, 0.0f);
	result.x = v1.y * v2.z - v1.z * v2.y;
	result.y = v1.z * v2.x - v1.x * v2.z;
	result.z = v1.x * v2.y - v1.y * v2.x;
	return result;
}