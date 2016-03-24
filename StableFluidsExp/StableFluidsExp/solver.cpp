#include "Helpers.h"

#define IX(i,j) ((i)*(N+2)+(j))
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
#define END_FOR }}

void add_source(int N, float * x, float * s, float dt)
{
	int i, size = (N + 2)*(N + 2);
	for (i = 0; i < size; i++) x[i] += dt*s[i];
}

void set_bnd(int N, int b, float * x)
{
	int i;

	for (i = 1; i <= N; i++) {
		x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
		x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
		x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
		x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
	}
	x[IX(0, 0)] = 0.5f*(x[IX(1, 0)] + x[IX(0, 1)]);
	x[IX(0, N + 1)] = 0.5f*(x[IX(1, N + 1)] + x[IX(0, N)]);
	x[IX(N + 1, 0)] = 0.5f*(x[IX(N, 0)] + x[IX(N + 1, 1)]);
	x[IX(N + 1, N + 1)] = 0.5f*(x[IX(N, N + 1)] + x[IX(N + 1, N)]);
}

void Gauss_Seidel(int N, int b, float * x, float * x0, float a, float c)
{
	int i, j, k;

	for (k = 0; k < 30; k++) {
		FOR_EACH_CELL
			x[IX(i, j)] = (x0[IX(i, j)] + a*(x[IX(i - 1, j)] + x[IX(i + 1, j)] + x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
		END_FOR
			set_bnd(N, b, x);
	}
}

void Jacobi_solve(int N, int b, float * x, float * x0, float a, float c)
{
	int i, j, k;
	int size = (N + 2) * (N + 2);
	float* aux = (float*)malloc(size*sizeof(float));

	for (k = 0; k < 50; k++)
	{
		FOR_EACH_CELL
			aux[IX(i, j)] = (x0[IX(i, j)] + a*(x[IX(i - 1, j)] + x[IX(i + 1, j)] + x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
		END_FOR

			FOR_EACH_CELL
			x[IX(i, j)] = aux[IX(i, j)];
		END_FOR
			set_bnd(N, b, x);
	}

	free(aux);
}

void diffuse(int N, int b, float * x, float * x0, float diff, float dt)
{
	float a = dt*diff*N*N;
	Jacobi_solve(N, b, x, x0, a, 1 + 4 * a);
}

float lerp(float t, float x0, float x1)
{
	return (1 - t) * x0 + t * x1;
}

void advect(int N, int b, float * d, float * d0, float * u, float * v, float dt)
{
	int i, j, i0, j0, i1, j1;
	float x, y, s1, t1, dt0;

	dt0 = dt*N;
	FOR_EACH_CELL
		x = i - dt0*u[IX(i, j)];
		y = j - dt0*v[IX(i, j)];

		if (x<0.5f) x = 0.5f;
		if (x>N + 0.5f) x = N + 0.5f;

		i0 = (int)x;
		i1 = i0 + 1;

		if (y<0.5f) y = 0.5f;
		if (y>N + 0.5f) y = N + 0.5f;

		j0 = (int)y;
		j1 = j0 + 1;

		s1 = x - i0;
		t1 = y - j0;

		float top_x_dir_lerp = lerp(s1, d0[IX(i0, j0)], d0[IX(i1, j0)]);
		float bottom_x_dir_lerp = lerp(s1, d0[IX(i0, j1)], d0[IX(i1, j1)]);
		d[IX(i, j)] = lerp(t1, top_x_dir_lerp, bottom_x_dir_lerp);
	END_FOR
		set_bnd(N, b, d);
}

void advect(int N, int b, float * d, float * d0, float * k, float * k0, float * u, float * v, float dt)
{
	int i, j, i0, j0, i1, j1;
	float x, y, s1, t1, dt0;

	dt0 = dt*N;
	FOR_EACH_CELL
		x = i - dt0*u[IX(i, j)];
		y = j - dt0*v[IX(i, j)];

		if (x<0.5f) x = 0.5f;
		if (x>N + 0.5f) x = N + 0.5f;

		i0 = (int)x;
		i1 = i0 + 1;

		if (y<0.5f) y = 0.5f;
		if (y>N + 0.5f) y = N + 0.5f;

		j0 = (int)y;
		j1 = j0 + 1;

		s1 = x - i0;
		t1 = y - j0;

		d[IX(i, j)] = lerp(s1, lerp(t1, d0[IX(i0, j0)], d0[IX(i0, j1)]), lerp(t1, d0[IX(i1, j0)], d0[IX(i1, j1)]));
		k[IX(i, j)] = lerp(s1, lerp(t1, k0[IX(i0, j0)], k0[IX(i0, j1)]), lerp(t1, k0[IX(i1, j0)], k0[IX(i1, j1)]));
	END_FOR

	set_bnd(N, b, d);
	set_bnd(N, b, k);
}

void project(int N, float * u, float * v, float * p, float * div)
{
	int i, j, iter = 0;

	computeDivergence_unifrom(N, u, v, div);
	set_bnd(N, 0, div);
	zeros(N, p);
	set_bnd(N, 0, p);

	Jacobi_solve(N, 0, p, div, 1, 4);

	FOR_EACH_CELL
		u[IX(i, j)] -= 0.5f*N*(p[IX(i + 1, j)] - p[IX(i - 1, j)]);
		v[IX(i, j)] -= 0.5f*N*(p[IX(i, j + 1)] - p[IX(i, j - 1)]);
	END_FOR
	set_bnd(N, 0, u); set_bnd(N, 0, v);
}

void add_force(int N, float dt, float * u, float * v, float * fx, float * fy)
{
	int i, j;
	FOR_EACH_CELL
		v[IX(i, j)] += dt * fy[IX(i, j)];
		u[IX(i, j)] += dt * fx[IX(i, j)];
	END_FOR
}

void dens_step(int N, float * x, float * x0, float * u, float * v, float diff, float dt)
{
	add_source(N, x, x0, dt);
	SWAP(x0, x); diffuse(N, 0, x, x0, diff, dt);
	SWAP(x0, x); advect(N, 0, x, x0, u, v, dt);
}

void vel_step(int N,
			  float * fx, float * fy,
			  float * psi,
			  float * du, float * dv,
			  float * wn, float *dw,
			  float * w_bar, float * w_star,
			  float * u, float * v, 
			  float * u0, float * v0, 
			  float visc, float dt)
{
	//This is time consuming but naive
	if (system("CLS")) system("clear");

	// IVOCK scheme
	zeros(N, wn);
	zeros(N, w_bar);
	zeros(N, w_star);
	zeros(N, dw);
	zeros(N, psi);
	zeros(N, du);
	zeros(N, dv);

	add_source(N, u, u0, dt); 
	add_source(N, v, v0, dt);
	
	SWAP(u0, u);
	SWAP(v0, v);
	diffuse(N, 0, u, u0, visc, dt);
	diffuse(N, 0, v, v0, visc, dt);
	project(N, u, v, u0, v0);
	SWAP(u0, u); 
	SWAP(v0, v);
	computeCurls_uniform(N, wn, u0, v0);
	set_bnd(N, 0, wn);
	//cout << endl << "Curl field from previous time step velocity field" << endl;
	//displayField(N + 2, N + 2, wn);
	advect(N, 0, w_bar, wn, u0, v0, dt);
	//cout << endl << "Curl field advected" << endl;
	//displayField(N + 2, N + 2, w_bar);
	advect(N, 0, u, u0, v, v0, u0, v0, dt);
	computeCurls_uniform(N, w_star, u, v);
	set_bnd(N, 0, w_star);
	//cout << endl << "Curl field from the advected velocity field" << endl;
	//displayField(N + 2, N + 2, w_star);
	linear_combine_sub(N, dw, w_bar, w_star);
	scaler(N, dw, -1.0f);
	set_bnd(N, 0, dw);
	//cout << endl << "Curl difference" << endl;
	//displayField(N + 2, N + 2, dw);
	Jacobi_solve(N, 0, psi, dw, 1, 4);
	//cout << endl << "Stream function (Psi)" << endl;
	//displayField(N + 2, N + 2, psi);
	find_vector_potential_2D(N, du, dv, psi);
	set_bnd(N, 0, du);
	set_bnd(N, 0, dv);
	//linear_combine_add(N, u, u, du);
	//linear_combine_add(N, v, v, dv);
	//set_bnd(N, 0, u);
	//set_bnd(N, 0, v);
	project(N, u, v, u0, v0);
}