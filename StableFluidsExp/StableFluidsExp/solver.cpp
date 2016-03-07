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

	for (k = 0; k < 20; k++) {
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

	for (k = 0; k < 20; k++)
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

void advect_beta(int N, int b, float * d, float * d0, float * k, float * k0, float * u, float * v, float dt)
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

		d[IX(i, j)] = lerp(t1, lerp(s1, d0[IX(i0, j0)], d0[IX(i1, j0)]), lerp(s1, d0[IX(i0, j1)], d0[IX(i1, j1)]));
		k[IX(i, j)] = lerp(t1, lerp(s1, k0[IX(i0, j0)], k0[IX(i1, j0)]), lerp(s1, k0[IX(i0, j1)], k0[IX(i1, j1)]));
	END_FOR

	set_bnd(N, b, d);
	set_bnd(N, b, k);
}

void project(int N, float * u, float * v, float * p, float * div)
{
	int i, j;

	computeDivergence_unifrom(N, u, v, div);
	FOR_EACH_CELL
		//div[IX(i, j)] = -0.5f*(u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]) / N;
		p[IX(i, j)] = 0;
	END_FOR
	set_bnd(N, 0, div); set_bnd(N, 0, p);

	cout << endl << "Divergence computed from velocity field" << endl;
	displayField(N + 2, N + 2, div);

	Gauss_Seidel(N, 0, p, div, 1, 4);

	FOR_EACH_CELL
		u[IX(i, j)] -= 0.5f*N*(p[IX(i + 1, j)] - p[IX(i - 1, j)]);
		v[IX(i, j)] -= 0.5f*N*(p[IX(i, j + 1)] - p[IX(i, j - 1)]);
	END_FOR
		set_bnd(N, 0, u); set_bnd(N, 0, v);
}

void dens_step(int N, float * x, float * x0, float * u, float * v, float diff, float dt)
{
	add_source(N, x, x0, dt);
	SWAP(x0, x); diffuse(N, 0, x, x0, diff, dt);
	SWAP(x0, x); advect(N, 0, x, x0, u, v, dt);
}

void vel_step(int N, float * w, float * w0, float * u, float * v, float * u0, float * v0, float visc, float dt)
{
	add_source(N, u, u0, dt); add_source(N, v, v0, dt);
	set_bnd(N, 0, u); set_bnd(N, 0, v);
	cout << endl << "Velocity field before stable fluids advction scheme" << endl;
	displayVectorField(N + 2, N + 2, u, v);

	cout << endl << "Curl field is computed from the velocity field" << endl;
	computeCurls_uniform(N + 2, w0, u, v);
	displayField(N + 2, N + 2, w0);

	SWAP(u0, u); SWAP(v0, v);
	diffuse(N, 0, u, u0, visc, dt);
	diffuse(N, 0, v, v0, visc, dt);
	advect_beta(N, 0, u, u0, v, v0, u0, v0, dt);
	cout << endl << "Velocity field after stable fluids advection scheme" << endl;
	displayVectorField(N + 2, N + 2, u, v);

	project(N, u, v, u0, v0);
	cout << endl << "Velocity field after pressure correction" << endl;
	displayVectorField(N + 2, N + 2, u, v);
	//SWAP(u0, u); SWAP(v0, v);
	//advect_beta(N, 1, u, u0, v, v0, u0, v0, dt);
	//cout << "Velocity field after stable fluids advection scheme" << endl;
	//displayVectorField(N + 2, N + 2, u, v);
	//project(N, u, v, u0, v0);
}

