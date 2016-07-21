#include <iostream>
#include "cptstd.hpp"
#include "matrix.hpp"
using namespace std;
using namespace cpt;

float accuracy = 0.001;        // desired relative accuracy in solution
int L = 4;                     // number of interior points in each dimension
int n_smooth = 5;               // number of pre and post smoothing iterations

Matrix<float, 2> psi(L + 2, L + 2), // solution to be found
psi_new(L + 2, L + 2),          // approximate solution after 1 iteration
rho(L + 2, L + 2);              // given source function

float h = 1.0 / (L + 1);         // step size
int steps;                      // number of iteration steps

void initialize()
{
	// check that L is a power of 2 as required by multigrid
	int power_of_2 = 1;
	while (power_of_2 < L)
		power_of_2 *= 2;
	if (power_of_2 != L) {
		L = power_of_2;
		cout << " Setting L = " << L << " (must be a power of 2)" << endl;
	}

	// create (L+2)x(L+2) matrices and zero them
	psi = psi_new = rho = Matrix<float, 2>(L + 2, L + 2);

	h = 1 / float(L + 1);      // assume physical size in x and y = 1
	float q = 10;              // point charge
	int i = L / 2;              // center of lattice
	rho[i][i] = q / (h * h);    // charge density

	for (int i = 1; i <= L; i++)
	{
		for (int j = 1; j <= L; j++)
			rho[i][j] = 0.0f;
	}

	rho[3][3] = 0.05f;
	rho[3][4] = 1.2f;

	steps = 0;
}

void Gauss_Seidel_I(Matrix<float, 2>& u, const Matrix<float, 2>& f, int iterations)
{
	// copy V to V_new
	// u = f;

	// Gauss-Seidel update in place
	for (int k = 0; k != iterations; k++)
	for (int i = 1; i <= L; i++)
	for (int j = 1; j <= L; j++)
		u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] +
		u[i][j - 1] + u[i][j + 1] +
		f[i][j]);
}

void Gauss_Seidel(float h, Matrix<float, 2>& u, const Matrix<float, 2>& f)
{
	int L = u.dim1() - 2;

	// use checkerboard updating
	for (int color = 0; color < 2; color++)
	for (int i = 1; i <= L; i++)
	for (int j = 1; j <= L; j++)
	if ((i + j) % 2 == color)
		u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] +
		u[i][j - 1] + u[i][j + 1] +
		f[i][j]);
}

void two_grid(float h, Matrix<float, 2>& u, Matrix<float, 2>& f)
{
	// solve exactly if there is only one interior point
	int L = u.dim1() - 2;
	if (L == 1) {
		u[1][1] = 0.25 * (u[0][1] + u[2][1] + u[1][0] + u[1][2] +
			f[1][1]);
		return;
	}

	// do a few pre-smoothing Gauss-Seidel steps
	for (int i = 0; i < n_smooth; i++)
		Gauss_Seidel(h, u, f);

	// find the residual
	Matrix<float, 2> r(L + 2, L + 2);
	for (int i = 1; i <= L; i++)
	for (int j = 1; j <= L; j++)
		r[i][j] = f[i][j] +
		(u[i + 1][j] + u[i - 1][j] +
		u[i][j + 1] + u[i][j - 1] - 4 * u[i][j]);

	// restrict residual to coarser grid
	int L2 = L / 2;
	Matrix<float, 2> R(L2 + 2, L2 + 2);
	for (int I = 1; I <= L2; I++) {
		int i = 2 * I - 1;
		for (int J = 1; J <= L2; J++) {
			int j = 2 * J - 1;
			R[I][J] = 0.25 * (r[i][j] + r[i + 1][j] + r[i][j + 1] +
				r[i + 1][j + 1]);
		}
	}

	// initialize correction V on coarse grid to zero
	Matrix<float, 2> V(L2 + 2, L2 + 2);

	// call twoGrid recursively
	float H = 2 * h;

	two_grid(H, V, R);

	// prolongate V to fine grid using simple injection
	Matrix<float, 2> v(L + 2, L + 2);
	for (int I = 1; I <= L2; I++) {
		int i = 2 * I - 1;
		for (int J = 1; J <= L2; J++) {
			int j = 2 * J - 1;
			v[i][j] = v[i + 1][j] = v[i][j + 1] = v[i + 1][j + 1] = V[I][J];
		}
	}

	// correct u
	for (int i = 1; i <= L; i++)
	for (int j = 1; j <= L; j++)
		u[i][j] += v[i][j];

	// do a few post-smoothing Gauss-Seidel steps
	for (int i = 0; i < n_smooth; i++)
		Gauss_Seidel(h, u, f);
}

float relative_error()
{
	float error = 0;           // average relative error per lattice point
	int n = 0;                  // number of non-zero differences

	for (int i = 1; i <= L; i++)
	for (int j = 1; j <= L; j++) {
		if (psi_new[i][j] != 0.0)
		if (psi_new[i][j] != psi[i][j]) {
			error += abs(1 - psi[i][j] / psi_new[i][j]);
			++n;
		}
	}
	if (n != 0)
		error /= n;

	return error;
}

int main()
{
	Matrix<float, 2> pressure(L + 2, L + 2);
	Matrix<float, 2> pressure_gs_results(L + 2, L + 2);
	Matrix<float, 2> div(L + 2, L + 2);
	float h = 1.0f / (L + 1);

	for (int i = 0; i != L + 2; i++){
		for (int j = 0; j != L + 2; j++){
			pressure[i][j] = div[i][j] = pressure_gs_results[i][j] = 0.0f;
		}
	}

	div[3][3] = 0.05f;
	div[3][4] = 1.2f;

	cout << std::setprecision(4) << setiosflags(ios::fixed);
	cout << "div field" << endl;
	for (int i = 0; i < L + 2; i++){
		for (int j = 0; j < L + 2; j++){
			cout << div[i][j] << ", ";
			if (j == L + 1)
				cout << endl;
		}
	}

	cout << "Multigrid solve" << endl;
	two_grid(h, pressure, div);

	cout << "GS solve" << endl;
	Gauss_Seidel_I(pressure_gs_results, div, 30);

	cout << "pressure field" << endl;
	for (int i = 0; i < L + 2; i++){
		for (int j = 0; j < L + 2; j++){
			cout << pressure[i][j] << ", ";
			if (j == L + 1)
				cout << endl;
		}
	}

	cout << "pressure field from GS solve" << endl;
	for (int i = 0; i < L + 2; i++){
		for (int j = 0; j < L + 2; j++){
			cout << pressure_gs_results[i][j] << ", ";
			if (j == L + 1)
				cout << endl;
		}
	}

	cout << "pressure solved" << endl;

	system("Pause");
	return 0;
}