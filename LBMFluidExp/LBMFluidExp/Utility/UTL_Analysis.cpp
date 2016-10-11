/*
* UTL_Analysis.cpp
* Xue Cui
* 14th 2016, Jan
* Desc:
*/
#include "UTL_Analysis.h"

namespace VFXEpoch
{
	namespace Analysis
	{
		// How to calculate curl in the uniform grid£¿
		// Localtion for storing the curl: Omega
		// 1.Calculate u component gradient from Cell6 & Cell4
		// 2.Calculate v component gradient from Cell8 & Cell2
		// 3.Subtract du from dv to get curl (dv - du).
		// Omega = dv - du;
		/*
			-------------------------
			|		| 		| 		|
			| Cell1	| Cell2 | Cell3 |
			|		|		|		|
			-------------------------
			|		|		|		|
			| Cell4	|  Curl | Cell6	|
			|		|		|		|
			-------------------------
			|		|		|		|
			| Cell7 | Cell8 | Cell9 |
			|		|		|		|
			-------------------------
		*/
		void computeCurl_uniform(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DVector2DfField ref) {
			float hx = dest.getDx();
			float hy = dest.getDy();
			for (int i = 1; i != dest.getDimY() - 1; i++){
				for (int j = 1; j != dest.getDimX() - 1; j++){
					float du, dv, curl;
					du = 0.5f * (ref.getData(i + 1, j).m_x - ref.getData(i - 1, j).m_x) / hx;
					dv = 0.5f * (ref.getData(i, j + 1).m_y - ref.getData(i, j - 1).m_y) / hy;
					curl = dv - du;
					dest.setData(curl, i, j);
				}
			}
		}

		void computeCurl_uniform(VFXEpoch::Grid3DVector3DfField& dest, VFXEpoch::Grid3DVector3DfField ref)
		{
			// TODO: Process 3D curl calculation in the uniform grid

		}
		
		// How to calculate curl in the mac grid£¿
		// Omega = Nabla x velocity field;
		// du = u2 - u1;
		// dv = v2 - v1;
		// Omega = du - dv;
		//
		/*
		-------------------------
		|			|			|
		|			u1
		|			|			|
		----v1-----Curl----v2----
		|			|			|
		|			u2			|
		|			|			|
		-------------------------
		*/
		void computeCurl_mac(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DfScalarField u, VFXEpoch::Grid2DfScalarField v) {
			for (int i = 1; i != dest.getDimY() - 1; i++){
				for (int j = 1; j != dest.getDimX() - 1; j++){
					float du, dv;
					float curl;
					du = (u.getData(i, j) - u.getData(i - 1, j)) / u.getDy();
					dv = (v.getData(i, j) - v.getData(i, j - 1)) / v.getDx();
					curl = dv - du;
					dest.setData(curl, i, j);
				}
			}
		}

		void computeCurl_mac(VFXEpoch::Grid3DVector3DfField& dest, VFXEpoch::Grid3DfScalarField u, VFXEpoch::Grid3DfScalarField v, VFXEpoch::Grid3DfScalarField w)
		{
			// TODO: Process 3D curl calculation in the uniform grid
		}

		void computeGradient_uniform(VFXEpoch::Grid2DVector2DfField& dest, VFXEpoch::Grid2DfScalarField ref) {
			// TODO: Compute gradients of a scalar field.
			for (int i = 1; i != dest.getDimY() - 1; i++){
				for (int j = 1; j != dest.getDimX() - 1; j++){
					VFXEpoch::Vector2Df data(ref.getData(i + 1, j) - ref.getData(i - 1, j), ref.getData(i, j + 1) - ref.getData(i, j - 1));
					dest.setData(data, i, j);
				}
			}
		}

		void
		computeGradient_mac(VFXEpoch::Grid2DVector2DfField& dest, VFXEpoch::Grid2DfScalarField ref)
		{
			// TODO: Compute gradients on staggered grid.
		}

		void
		find_vector_from_vector_potential_2D(VFXEpoch::Grid2DVector2DfField& u, VFXEpoch::Grid2DfScalarField psi)
		{
			// TODO: Calculate vector potential from the given function (Psi for stream function)
			float dpsi_dy, dpsi_dx;
			int Nx = u.getDimX() - 2;
			int Ny = u.getDimY() - 2;
			u.zeroVectors();
						
			for (int i = 1; i != u.getDimY() - 1; i++){
				for (int j = 1; j != u.getDimX() - 1; j++){
					dpsi_dy = 0.5f * (psi.getData(i + 1, j) - psi.getData(i - 1, j)) * Ny;
					dpsi_dx = 0.5f * (psi.getData(i, j + 1) - psi.getData(i, j - 1)) * Nx;
					u.dx = dpsi_dy;
					u.dy = -dpsi_dx;
				}
			}
		}

		void computeDivergence_uniform(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DVector2DfField ref, float scale) {
			for (int i = 1; i != dest.getDimY() - 1; i++){
				for (int j = 1; j != dest.getDimX() - 1; j++){
					float div = 0.0f;
					float du, dv;

					int Nx = ref.getDimX() - 2;
					int Ny = ref.getDimY() - 2;

					du = scale * 0.5f * (ref.getData(i + 1, j).m_x - ref.getData(i - 1, j).m_x) / Nx;
					dv = scale * 0.5f * (ref.getData(i, j + 1).m_y - ref.getData(i, j - 1).m_y) / Ny;
					div = du + dv;
					dest.setData(div, i, j);
				}
			}
		}

		void computeDivergence_mac(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DfScalarField u, VFXEpoch::Grid2DfScalarField v) {
			// TODO: Compute divergence on staggered grid
		}
	}
}