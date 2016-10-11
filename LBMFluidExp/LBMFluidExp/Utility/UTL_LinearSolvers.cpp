/*
 * UTL_LinearSolvers.cpp
 * Xue Cui
 * 8th Mar, 2016
 * Desc:
 *
*/
#include "UTL_LinearSolvers.h"

namespace VFXEpoch
{
	namespace LinearSolver
	{
		void 
		GSSolve(VFXEpoch::Grid2DfScalarField& x, VFXEpoch::Grid2DfScalarField x0, VFXEpoch::BndConditionPerEdge b[], float coefMatrixAElement, float c, int iterations)
		{
			for (int m = 0; m != iterations; m++)
			{
				for (int i = 1; i != x.getDimY() - 1; i++)
				{
					for (int j = 1; j != x.getDimX() - 1; j++)
					{
						float sumFluxout = x.getData(i + 1, j) + x.getData(i - 1, j) + x.getData(i, j + 1) + x.getData(i, j - 1);
						sumFluxout *= coefMatrixAElement;
						float curApproximation = x0.getData(i, j) + sumFluxout;
						curApproximation /= c;
						x.setData(curApproximation, i, j);
					}
				}
				// TODO: Boundary values need to be set
				x.setBoundaries(b[0].boundaryType, b[0].side);
				x.setBoundaries(b[1].boundaryType, b[1].side);
				x.setBoundaries(b[2].boundaryType, b[2].side);
				x.setBoundaries(b[3].boundaryType, b[3].side);
				x.setBoundariesOnCorners();
			}
		}

		void 
		JacobiSolve(VFXEpoch::Grid2DfScalarField& x, VFXEpoch::Grid2DfScalarField x0, VFXEpoch::BndConditionPerEdge b[], float coefMatrixAElement, float c, int iterations)
		{
			// Extra memory to avoid covering previous value
			VFXEpoch::Grid2DfScalarField auxiliary(x0.getDimY(), x0.getDimX(), (float)(1 / x0.getDimX()), (float)(1 / x0.getDimY()));

			for (int m = 0; m != iterations; m++)
			{
				for (int i = 1; i != x.getDimY() - 1; i++)
				{
					for (int j = 1; j != x.getDimX() - 1; j++)
					{
						float sumFluxout = x.getData(i + 1, j) + x.getData(i - 1, j) + x.getData(i, j + 1) + x.getData(i, j - 1);
						sumFluxout *= coefMatrixAElement;
						float curApproximation = x0.getData(i, j) + sumFluxout;
						curApproximation /= c;
						auxiliary.setData(curApproximation, i, j);
					}
				}

				for (int i = 1; i != x.getDimY() - 1; i++)
				{
					for (int j = 1; j != x.getDimX() - 1; j++)
					{
						x.setData(auxiliary.getData(i, j), i, j);
					}
				}

				// TODO: Boundary values need to be set
				x.setBoundaries(b[0].boundaryType, b[0].side);
				x.setBoundaries(b[1].boundaryType, b[1].side);
				x.setBoundaries(b[2].boundaryType, b[2].side);
				x.setBoundaries(b[3].boundaryType, b[3].side);
				x.setBoundariesOnCorners();
			}
		}
	}
}