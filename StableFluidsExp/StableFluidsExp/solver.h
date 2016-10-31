#ifndef _SOLVER_H_
#define _SOLVER_H_
#include "helpers.h"

void
MoveScalarProperties(int N, float * x, float * x0, float * u, float * v, float diff, float dt);

void 
IVOCKAdvance(int N,
			 Particle* particles, int num_particles,
			 float * fx, float * fy,
			 float * psi, float * du, float * dv, float * wn, float *dw, float * w_bar, float * w_star,
			 float * u, float * v, float * u0, float * v0,
			 float * t, float * t0,
			 float visc,
			 float dt);

void
LBMAdvance(float * f0,
		  float * f1, float * f2, float * f3, float * f4,
		  float * f5, float * f6, float * f7, float * f8,
		  int N, float tau, float * out_u, float * out_v,
		  Particle* particles, int num_particles, float dt);

void
project(int N, float * u, float * v, float * p, float * div);

void
set_boundaries(int N, int b, float * x);

void
add_gravity(int N, float dt, float grav, float * field);

#endif _SOLVER_H_