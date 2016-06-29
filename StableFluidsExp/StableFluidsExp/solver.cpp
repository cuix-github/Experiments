#include "Helpers.h"

#define IX(i,j) ((i)*(N+2)+(j))
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
#define END_FOR }}

void
add_source(int N, float * x, float * s, float dt){
	int i, size = (N + 2)*(N + 2);
	for (i = 0; i < size; i++) x[i] += dt*s[i];
}

void
set_boundaries(int N, int b, float * x){
	for (int i = 1; i <= N; i++) {
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

void
Gauss_Seidel_solve(int N, int b, float * x, float * x0, float a, float c, int iterations){
	int i, j, k;

	double h = 1.0f / N;

	for (k = 0; k < iterations; k++) {
		FOR_EACH_CELL
			x[IX(i, j)] = (x0[IX(i, j)] + a *
			(x[IX(i - 1, j)] + x[IX(i + 1, j)] +
			x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
		END_FOR
		set_boundaries(N, b, x);
	}
}

void
Jacobi_solve(int N, int b, float * x, float * x0, float a, float c, int iterations){
	int i, j, k;
	int size = (N + 2) * (N + 2);
	float* aux = (float*)malloc(size*sizeof(float));
	double h = 1.0f / (N);

	for (k = 0; k < iterations; k++)
	{
		FOR_EACH_CELL
			aux[IX(i, j)] = (x0[IX(i, j)] * 1.0f / (h * h) + a *
			(x[IX(i - 1, j)] + x[IX(i + 1, j)] +
			x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
		END_FOR

			FOR_EACH_CELL
			x[IX(i, j)] = aux[IX(i, j)];
		END_FOR
		set_boundaries(N, b, x);
	}

	free(aux);
}

void
diffuse(int N, int b, float * x, float * x0, float diff, float dt){
	float a = dt*diff*N*N;
	Gauss_Seidel_solve(N, b, x, x0, a, 1 + 4 * a, 30);
}

vec2
get_velocity(int N, const vec2& position, float * u, float * v){
	float u0, v0;
	u0 = interpolate(N, position.x * N - 0.0f, position.y * N - 0.5f, u);
	v0 = interpolate(N, position.x * N - 0.5f, position.y * N - 0.0f, v);
	return vec2(u0, v0);
}

// Mark:
// TODO: Problems make the particles move to wrong direction
// Desc:
// Runge-Kutta 2nd order integration for ODEs
vec2 rk2(int N, float * u, float * v, const vec2& position, float dt){
	vec2 vel = get_velocity(N, position, u, v);
	vel = get_velocity(N, vec2(position.x + 0.5 * dt * vel.x, position.y + 0.5 * dt * vel.y), u, v);
	return vec2(position.x + dt * vel.x, position.y + dt * vel.y);
}

void
particles_advector(int N, float * u, float * v, Particle* particles, int num_particles, float dt){
	vec2 pos(0.0f, 0.0f);
	for (int i = 0; i != num_particles; i++){
		pos.x = pos.y = 0.0f;
		pos = rk2(N, u, v, vec2(particles[i].x, particles[i].y), dt);
		particles[i].x = pos.x; particles[i].y = pos.y;
		particles[i].vel = get_velocity(N, vec2(particles[i].x, particles[i].y), u, v);
	}
}

// Problems: Blow up : (
// 2nd order Runge-Kutta ODEs integrator for advection bugs
void
vector_advector_rk2(int N, float * u, float * u0, float * v, float * v0, float * u_tmp, float * v_tmp, float dt){
	// TODO: Runge Kutta integrator (2nd order)
}

void
scalar_advector(int N, float * d, float * d0, float * u, float * v, float dt){
	int i, j, i0, j0, i1, j1;
	float x, y, s1, t1, dt0;

	dt0 = dt*N;
	FOR_EACH_CELL
		x = i - dt0*u[IX(i, j)];
		y = j - dt0*v[IX(i, j)];

		if (x < 0.5f) x = 0.5f;
		if (x > N + 0.5f) x = N + 0.5f;

		i0 = (int)x;
		i1 = i0 + 1;

		if (y < 0.5f) y = 0.5f;
		if (y > N + 0.5f) y = N + 0.5f;

		j0 = (int)y;
		j1 = j0 + 1;

		s1 = x - i0;
		t1 = y - j0;

		float top_x_dir_lerp = lerp(s1, d0[IX(i0, j0)], d0[IX(i1, j0)]);
		float bottom_x_dir_lerp = lerp(s1, d0[IX(i0, j1)], d0[IX(i1, j1)]);
		d[IX(i, j)] = lerp(t1, top_x_dir_lerp, bottom_x_dir_lerp);
	END_FOR
}

void
vector_advector(int N, float * d, float * d0, float * k, float * k0, float * u, float * v, float dt){
	int i, j, i0, j0, i1, j1;
	float x, y, s1, t1, dt0;

	// TODO: Modify this piece of code to archive 2nd Order RK
	dt0 = dt*N;
	FOR_EACH_CELL
		x = i - dt0*u[IX(i, j)];
		y = j - dt0*v[IX(i, j)];

		if (x < 0.5f) x = 0.5f;
		if (x > N + 0.5f) x = N + 0.5f;

		i0 = (int)x;
		i1 = i0 + 1;

		if (y < 0.5f) y = 0.5f;
		if (y > N + 0.5f) y = N + 0.5f;

		j0 = (int)y;
		j1 = j0 + 1;

		s1 = x - i0;
		t1 = y - j0;

		d[IX(i, j)] = lerp(s1,
			lerp(t1, d0[IX(i0, j0)], d0[IX(i0, j1)]),
			lerp(t1, d0[IX(i1, j0)], d0[IX(i1, j1)]));

		k[IX(i, j)] = lerp(s1,
			lerp(t1, k0[IX(i0, j0)], k0[IX(i0, j1)]),
			lerp(t1, k0[IX(i1, j0)], k0[IX(i1, j1)]));
	END_FOR
}

void
project(int N, float * u, float * v, float * p, float * div){
	int i, j, iter = 0;
	zeros(N, p);

	computeDivergence_unifrom(N, u, v, div);
	scaler(N, div, -1.0f);
	Gauss_Seidel_solve(N, 0, p, div, 1, 4, 50);
	FOR_EACH_CELL
		u[IX(i, j)] -= 0.5f*N*(p[IX(i + 1, j)] - p[IX(i - 1, j)]);
		v[IX(i, j)] -= 0.5f*N*(p[IX(i, j + 1)] - p[IX(i, j - 1)]);
	END_FOR
	set_boundaries(N, 0, u);
	set_boundaries(N, 0, v);
}

void
add_force(int N, float dt, float * u, float * v, float * fx, float * fy){
	int i, j;
	FOR_EACH_CELL
		v[IX(i, j)] += dt * fy[IX(i, j)];
	u[IX(i, j)] += dt * fx[IX(i, j)];
	END_FOR
}

void
MoveScalarProperties(int N, float * x, float * x0, float * u, float * v, float diff, float dt){
	add_source(N, x, x0, dt);
	SWAP(x0, x); diffuse(N, 0, x, x0, diff, dt);
	SWAP(x0, x); scalar_advector(N, x, x0, u, v, dt);
}

void IVOCKAdvance(int N,
	Particle* particles, int num_particles,
	float * fx, float * fy,
	float * psi, float * du, float * dv, float * wn, float *dw, float * w_bar, float * w_star,
	float * u, float * v, float * u0, float * v0,
	float * t, float * t0,
	float visc,
	float dt){

	// Only for debug

	//// IVOCK advection
	zeros(N, wn);
	zeros(N, w_bar);
	zeros(N, w_star);
	zeros(N, dw);
	zeros(N, psi);
	zeros(N, du);
	zeros(N, dv);

	// Gravity
	int size = (N + 2) * (N + 2);
	float *g = (float*)malloc(size*sizeof(float));
	zeros(N, g);
	for (int i = 1; i <= N; i++){
		for (int j = 1; j <= N; j++){
			g[IX(i, j)] = -0.098f;
		}
	}

	add_source(N, u, u0, dt);
	add_source(N, v, v0, dt);
	add_source(N, v, g, dt);
	add_source(N, t, t0, dt);

	particles_advector(N, u, v, particles, num_particles, dt);

	SWAP(u0, u);
	SWAP(v0, v);
	diffuse(N, 0, u, u0, visc, dt);
	diffuse(N, 0, v, v0, visc, dt);
	project(N, u, v, u0, v0);
	zeros(N, u0);
	zeros(N, v0);
	SWAP(u0, u);
	SWAP(v0, v);

	computeCurls_uniform(N, wn, u0, v0);
	scalar_advector(N, w_bar, wn, u0, v0, dt);
	vector_advector(N, u, u0, v, v0, u0, v0, dt);
	computeCurls_uniform(N, w_star, u, v);
	linear_combine_sub(N, dw, w_bar, w_star);
	scaler(N, dw, -1.0f);
	Jacobi_solve(N, 0, psi, dw, -1, -4, 50);
	find_vector_potential_2D(N, du, dv, psi);
	//linear_combine_add(N, u, u, du);
	//linear_combine_add(N, v, v, dv);
	free(g);
}

// Poisson Equation Laplace(Psi) = f(x);
// Ex.1
// 5x5 with interior field 3x3 and outter boundaries
/*
A =
-4	1	0	1	0	0	0	0	0
1	-4	1	0	1	0	0	0	0
0	1	-4	0	0	1	0	0	0
1	0	0	-4	1	0	1	0	0
0	1	0	1	-4	1	0	1	0
0	0	1	0	1	-4	0	0	1
0	0	0	1	0	0	-4	1	0
0	0	0	0	1	0	1	-4	1
0	0	0	0	0	1	0	1	-4

b' = 5 0 0 0 6 0 0 0 0

10 Times Gauss Seidel solution is:
Step 1					Step 2				Step 3				Step 4				Step 5
-1.25000000000000	-1.40625000000000	-1.64257812500000	-1.75708007812500	-1.81387329101563
-0.312500000000000	-0.785156250000000	-1.01416015625000	-1.12774658203125	-1.18442535400391
-0.0781250000000000	-0.304687500000000	-0.417968750000000	-0.474609375000000	-0.502929687500000
-0.312500000000000	-0.785156250000000	-1.01416015625000	-1.12774658203125	-1.18442535400391
-1.65625000000000	-2.10937500000000	-2.33593750000000	-2.44921875000000	-2.50585937500000
-0.433593750000000	-0.657714843750000	-0.770690917968750	-0.827293395996094	-0.855608940124512
-0.0781250000000000	-0.304687500000000	-0.417968750000000	-0.474609375000000	-0.502929687500000
-0.433593750000000	-0.657714843750000	-0.770690917968750	-0.827293395996094	-0.855608940124512
-0.216796875000000	-0.328857421875000	-0.385345458984375	-0.413646697998047	-0.427804470062256
The exact solution for this system is:
(-1.8705   -1.2411   -0.5313   -1.2411   -2.5625   -0.8839   -0.5313   -0.8839   -0.4420)

Step 6					Step 7					Step 8				Step 9				Step 10
-1.84221267700195	-1.85637521743774	-1.86345559358597	-1.86699566990137	-1.86876569408923
-1.21275043487549	-1.22691118717194	-1.23399133980274	-1.23753138817847	-1.23930140887387
-0.517089843750000	-0.524169921875000	-0.527709960937500	-0.529479980468750	-0.530364990234375
-1.21275043487549	-1.22691118717194	-1.23399133980274	-1.23753138817847	-1.23930140887387
-2.53417968750000	-2.54833984375000	-2.55541992187500	-2.55895996093750	-2.56072998046875
-0.869768500328064	-0.876848503947258	-0.880388533696532	-0.882158552063629	-0.883043561683735
-0.517089843750000	-0.524169921875000	-0.527709960937500	-0.529479980468750	-0.530364990234375
-0.869768500328064	-0.876848503947258	-0.880388533696532	-0.882158552063629	-0.883043561683735
-0.434884250164032	-0.438424251973629	-0.440194266848266	-0.441079276031815	-0.441521780841867
The exact solution for this system is:
(-1.8705   -1.2411   -0.5313   -1.2411   -2.5625   -0.8839   -0.5313   -0.8839   -0.4420)
*/