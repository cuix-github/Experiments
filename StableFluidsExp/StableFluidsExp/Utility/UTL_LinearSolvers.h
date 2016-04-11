/*
 * UTL_LinearSolvers.h
 * Xue Cui
 * 8th Mar, 2016
 * Desc:
 *
*/
#ifndef _UTL_LINEAR_SOLVERS_H_
#define _UTL_LINEAR_SOLVERS_H_
#include <math.h>
#include <iostream>
#include "UTL_Grid.h"

using namespace std;

namespace VFXEpoch
{
	namespace LinearSolver
	{
		void GSSolve(VFXEpoch::Grid2DfScalarField& x, VFXEpoch::Grid2DfScalarField x0, VFXEpoch::BndConditionPerEdge b[], float coefMatrixAElement, float c, int iterations);
		void JacobiSolve(VFXEpoch::Grid2DfScalarField& x, VFXEpoch::Grid2DfScalarField x0, VFXEpoch::BndConditionPerEdge b[], float coefMatrixAElement, float c, int iterations);
	}
}

#endif