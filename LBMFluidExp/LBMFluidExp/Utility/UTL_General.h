/*
* UTL_Matrix.h
* Xue Cui
* 13th Nov, 2015
* Desc:
*
*/
#ifndef _UTL_GENERAL_H_
#define _UTL_GENERAL_H_

#include <stdlib.h>
#include <iostream>
#include <random>
#include "UTL_Grid.h"

using namespace std;

namespace VFXEpoch
{
	enum class VECTOR_COMPONENTS
	{
		X, Y, Z
	};

	// Helpers
	float RandomF(float min, float max);
	int RandomI(int min, int max);
	float Lerp(float t, float x0, float x1);
	float Bilerp(float t, float s, float x0, float x1, float y0, float y1);
	void ExtractComponents(VFXEpoch::Grid2DfScalarField& component, VFXEpoch::Grid2DVector2DfField vectorField, VECTOR_COMPONENTS axis);
	void ExtractComponents(VFXEpoch::Grid2DfScalarField& component, VFXEpoch::Grid3DVector3DfField vectorField, VECTOR_COMPONENTS axis);
	void InsertComponents(VFXEpoch::Grid2DfScalarField component, VFXEpoch::Grid2DVector2DfField& vectorField, VECTOR_COMPONENTS axis);
	void InsertComponents(VFXEpoch::Grid2DfScalarField component, VFXEpoch::Grid3DVector3DfField& vectorField, VECTOR_COMPONENTS axis);
	void Zeros(VFXEpoch::Grid2DfScalarField& field);
	void Zeros(VFXEpoch::Grid2DVector2DfField& field);
	
	// Template Helpers
	template <class T>
	void Swap(T& a, T& b)
	{
		T c(a);
		a = b;
		b = c;
	}
}

#endif 