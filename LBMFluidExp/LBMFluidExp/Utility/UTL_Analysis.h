/*
* UTL_Analysis.h
* Xue Cui
* 2nd Jan, 2016
* Desc:
* Some general math operations
*/
#ifndef _UTL_ANALYSIS_H_
#define _UTL_ANALYSIS_H_

#include "UTL_Grid.h"

namespace VFXEpoch
{
	namespace Analysis
	{
		void 
		computeCurl_uniform(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DVector2DfField ref);
		
		void 
		computeCurl_uniform(VFXEpoch::Grid3DVector3DfField& dest, VFXEpoch::Grid3DVector3DfField ref);
		
		void 
		computeCurl_mac(VFXEpoch::Grid2DfScalarField& dest, 
						VFXEpoch::Grid2DfScalarField u, VFXEpoch::Grid2DfScalarField v);
		
		void 
		computeCurl_mac(VFXEpoch::Grid3DVector3DfField& dest, 
						VFXEpoch::Grid2DfScalarField u, VFXEpoch::Grid2DfScalarField v, VFXEpoch::Grid2DfScalarField w);
		
		void 
		computeGradient_uniform(VFXEpoch::Grid2DVector2DfField& dest, VFXEpoch::Grid2DfScalarField ref);
		
		void
		computeGradient_mac(VFXEpoch::Grid2DVector2DfField& dest, VFXEpoch::Grid2DfScalarField ref);

		void 
		computeDivergence_uniform(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DVector2DfField ref, float scale);

		void
		computeDivergence_mac(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DfScalarField u, VFXEpoch::Grid2DfScalarField v);
	}
}

#endif