#ifndef _SOLVER_H_
#define _SOLVER_H_
#include "helpers.h"

void
MoveDens(int N, float * x, float * x0, float * u, float * v, float diff, float dt);

void IVOCKAdvance(int N,
	Particle* particles, int num_particles,
	float * fx, float * fy,
	float * psi, float * du, float * dv, float * wn, float *dw, float * w_bar, float * w_star,
	float * u, float * v, float * u0, float * v0,
	float visc,
	float dt);

#endif _SOLVER_H_