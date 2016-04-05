#include "Helpers.h"

#define IX(i,j) ((i)*(N+2)+(j))
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
#define END_FOR }}

void add_source(int N, 
				float * x, 
				float * s, 
				float dt)
{
	int i, size = (N + 2)*(N + 2);
	for (i = 0; i < size; i++) x[i] += dt*s[i];
}

void 
set_bnd(int N, 
		int b, 
		float * x)
{
	for (int i = 1; i <= N; i++) {
		x[IX(0, i)]		= b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
		x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
		x[IX(i, 0)]		= b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
		x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
	}
	x[IX(0, 0)] = 0.5f*(x[IX(1, 0)] + x[IX(0, 1)]);
	x[IX(0, N + 1)] = 0.5f*(x[IX(1, N + 1)] + x[IX(0, N)]);
	x[IX(N + 1, 0)] = 0.5f*(x[IX(N, 0)] + x[IX(N + 1, 1)]);
	x[IX(N + 1, N + 1)] = 0.5f*(x[IX(N, N + 1)] + x[IX(N + 1, N)]);
}

void 
Gauss_Seidel(int N, 
			 int b, 
			 float * x, 
			 float * x0, 
			 float a, 
			 float c)
{
	int i, j, k;

	double h = 1 / double(N + 1);

	for (k = 0; k < 20; k++) {
		FOR_EACH_CELL
			x[IX(i, j)] = (x0[IX(i, j)] + a * 
						  ( x[IX(i - 1, j)] + x[IX(i + 1, j)] + 
						    x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
		END_FOR
		set_bnd(N, b, x);
		//cout << endl << "Step: " << k << endl;
		//displayField(N + 2, N + 2, x);
	}
}

void
Gauss_Seidel_Streamfunction(int N,
							int b,
							float * x,
							float * x0,
							float a,
							float c)
{
	int i, j, k;

	double h = 1 / double(N + 1);

	for (k = 0; k < 100; k++) {
		FOR_EACH_CELL
			x[IX(i, j)] = (x0[IX(i, j)] * h * h * 2.3f + a *
			(x[IX(i - 1, j)] + x[IX(i + 1, j)] +
			x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
		END_FOR
			set_bnd(N, b, x);
		//cout << endl << "Step: " << k << endl;
		//displayField(N + 2, N + 2, x);
	}
}

void 
Jacobi_solve(int N, 
			 int b, 
			 float * x, 
			 float * x0, 
			 float a, 
			 float c)
{
	int i, j, k;
	int size = (N + 2) * (N + 2);
	float* aux = (float*)malloc(size*sizeof(float));
	double h = 1 / (double)(N + 1);

	for (k = 0; k < 30; k++)
	{
		FOR_EACH_CELL
			aux[IX(i, j)] = (x0[IX(i, j)] + a * 
							( x[IX(i - 1, j)] + x[IX(i + 1, j)] + 
							  x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
		END_FOR

		FOR_EACH_CELL
			x[IX(i, j)] = aux[IX(i, j)];
		END_FOR
			set_bnd(N, b, x);

		//cout << endl << "Step: " << k << endl;
		//displayField(N + 2, N + 2, x);
	}

	free(aux);
}

void 
diffuse(int N, 
		int b, 
		float * x, 
		float * x0, 
		float diff, 
		float dt)
{
	float a = dt*diff*N*N;
	Jacobi_solve(N, b, x, x0, a, 1 + 4 * a);
}


void
get_barycentric(float x,
				int& i,
				float& f,
				int i_low,
				int i_high)
{
	float s = std::floor(x);
	i = (int)s;
	if (i<i_low){
		i = i_low;
		f = 0;
	}
	else if (i>i_high - 2){
		i = i_high - 2;
		f = 1;
	}
	else
		f = (float)(x - s);
}

float
interpolate(int N,
			float x,
			float y,
			float *field){
	int i, j;
	float fx, fy;
	get_barycentric(x, i, fx, 0, N);
	get_barycentric(y, j, fy, 0, N);
	return bilerp(fx, fy, field[IX(i, j)], field[IX(i, j + 1)], field[IX(i + 1, j)], field[IX(i + 1, j + 1)]);
}

void
getVel(int N,
	   float u_out,
	   float v_out,
	   float x,
	   float y,
	   float* u,
	   float* v){
	u_out = interpolate(N, x * N, y * N - 0.5f, u);
	v_out = interpolate(N, x * N - 0.5f, y * N, v);
}

void
trace_rk2(
		  int N,
		  float x_out, 
		  float y_out, 
		  float x, 
		  float y,
		  float * u,
		  float * v,
		  float dt){
	float _u, _v;
	getVel(N, _u, _v, x, y, u, v);
	getVel(N, _u, _v, x + 0.5f * dt * _u, y + 0.5f * dt * _v, u, v);
	x_out = x + 0.5f * dt * _u;
	y_out = y + 0.5f * dt * _v;
}

void
advect_particle(Particle* particles,
				int N,
				int numParticles,
				float * u, 
				float * v,
				float dt)
{
	for (int i = 0; i != numParticles; i++){
		trace_rk2(N, particles[i].x, particles[i].y, particles[i].x, particles[i].y, u, v, dt);
	}
}

void
advect_rk2(int N,
		   float dt,
		   float * u,
		   float * v,
		   float * u_temp,
		   float * v_temp
		   )
{
	float h = 1.0 / N;
	for (int i = 1; i <= N; i++){
		for (int j = 1; j <= N; j++){
			float x = i * h;
			float y = (j + 0.5) * h;
			trace_rk2(N, x, y, x, y, u, v, -dt);
			getVel(N, u_temp[IX(i, j)], 0, x, y, u, v);
		}
	}

	for (int i = 1; i <= N; i++){
		for (int j = 1; j <= N; j++){
			float x = (i + 0.5f) * h;
			float y = j * h;
			trace_rk2(N, x, y, x, y, u, v, -dt);
			getVel(N, 0, v_temp[IX(i, j)], x, y, u, v);
		}
	}

	u = u_temp;
	v = v_temp;
}

void 
advect(int N, 
	   int b, 
	   float * d, 
	   float * d0, 
	   float * u, 
	   float * v, 
	   float dt)
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

void 
advect(int N, 
	   int b, 
	   float * d, 
	   float * d0, 
	   float * k, 
	   float * k0, 
	   float * u, 
	   float * v, 
	   float dt)
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

		d[IX(i, j)] = lerp(s1, 
					  lerp(t1, d0[IX(i0, j0)], d0[IX(i0, j1)]), 
					  lerp(t1, d0[IX(i1, j0)], d0[IX(i1, j1)]));

		k[IX(i, j)] = lerp(s1, 
					  lerp(t1, k0[IX(i0, j0)], k0[IX(i0, j1)]), 
					  lerp(t1, k0[IX(i1, j0)], k0[IX(i1, j1)]));
	END_FOR

	set_bnd(N, b, d);
	set_bnd(N, b, k);
}

void 
project(int N, 
		float * u, 
		float * v, 
		float * p, 
		float * div)
{
	int i, j, iter = 0;

	computeDivergence_unifrom(N, u, v, div);
	zeros(N, p);
	Jacobi_solve(N, 0, p, div, 1, 4);
	FOR_EACH_CELL
		u[IX(i, j)] -= 0.5f*N*(p[IX(i + 1, j)] - p[IX(i - 1, j)]);
		v[IX(i, j)] -= 0.5f*N*(p[IX(i, j + 1)] - p[IX(i, j - 1)]);
	END_FOR
	set_bnd(N, 0, u); 
	set_bnd(N, 0, v);
}

void 
add_force(int N, 
		  float dt, 
		  float * u, 
		  float * v, 
		  float * fx, 
		  float * fy)
{
	int i, j;
	FOR_EACH_CELL
		v[IX(i, j)] += dt * fy[IX(i, j)];
		u[IX(i, j)] += dt * fx[IX(i, j)];
	END_FOR
}

void 
dens_step(int N, 
		  float * x, 
		  float * x0, 
		  float * u, 
		  float * v, 
		  float diff, 
		  float dt)
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
	//
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
	//set_bnd(N, 0, wn);
	//cout << endl << "Curl field from previous time step velocity field" << endl;
	//displayField(N + 2, N + 2, wn);
	advect(N, 0, w_bar, wn, u0, v0, dt);
	//cout << endl << "Curl field advected" << endl;
	//displayField(N + 2, N + 2, w_bar);
	advect(N, 0, u, u0, v, v0, u0, v0, dt);
	computeCurls_uniform(N, w_star, u, v);
	//set_bnd(N, 0, w_star);
	//cout << endl << "Curl field from the advected velocity field" << endl;
	//displayField(N + 2, N + 2, w_star);
	linear_combine_sub(N, dw, w_bar, w_star);
	scaler(N, dw, -1.0f);
	//set_bnd(N, 0, dw);
	//cout << endl << "Curl difference" << endl;
	//displayField(N + 2, N + 2, dw);
	//if (system("CLS")) system("clear");
	Gauss_Seidel_Streamfunction(N, 0, psi, dw, -1, -4);
	//cout << endl << "Stream function (Psi)" << endl;
	//displayField(N + 2, N + 2, psi);
	find_vector_potential_2D(N, du, dv, psi);
	//set_bnd(N, 0, du);
	//set_bnd(N, 0, dv);
	linear_combine_add(N, u, u, du);
	linear_combine_add(N, v, v, dv);
	//set_bnd(N, 0, u);
	//set_bnd(N, 0, v);
	project(N, u, v, u0, v0);
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