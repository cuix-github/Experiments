#include "../cpt/cptstd.hpp"
#include "../cpt/matrix.hpp"
#include "Helpers.h"
using namespace cpt;

int L = 50;                     // number of interior points in x and y
Matrix<double, 2> V(L + 2, L + 2),   // potential to be found
rho(L + 2, L + 2),              // given charge density
V_new(L + 2, L + 2);            // new potential after each step

double h;                       // lattice spacing
int steps;                      // number of iteration steps
double accuracy;                // desired accuracy in solution
double omega;                   // overrelaxation parameter

void initialize()
{
	int N = L + 2;
	V = rho = V_new = Matrix<double, 2>(N, N);

	h = 1 / double(L + 1);      // assume physical size in x and y = 1
	double q = 10;              // point charge
	int i = N / 2;              // center of lattice
	rho[i][i] = q / (h * h);    // charge density
	steps = 0;
}

void Jacobi() {
	// Jacobi algorithm for a single iterative step
	for (int i = 1; i <= L; i++)
	for (int j = 1; j <= L; j++)
		V_new[i][j] = 0.25 * (V[i - 1][j] + V[i + 1][j] +
		V[i][j - 1] + V[i][j + 1] +
		h * h * rho[i][j]);
}
double relative_error()
{
	double error = 0;           // average relative error per lattice point
	int n = 0;                  // number of non-zero differences

	for (int i = 1; i <= L; i++)
	for (int j = 1; j <= L; j++) {
		if (V_new[i][j] != 0)
		if (V_new[i][j] != V[i][j]) {
			error += abs(1 - V[i][j] / V_new[i][j]);
			++n;
		}
	}
	if (n != 0)
		error /= n;
	return error;
}

void Gauss_Seidel()
{
	// copy V to V_new
	V_new = V;

	// Gauss-Seidel update in place
	for (int i = 1; i <= L; i++)
	for (int j = 1; j <= L; j++)
		V_new[i][j] = 0.25 * (V_new[i - 1][j] + V_new[i + 1][j] +
		V_new[i][j - 1] + V_new[i][j + 1] +
		h * h * rho[i][j]);
}

void successive_over_relaxation()   // using red-black checkerboard updating
{
	// update even sites first
	for (int i = 1; i <= L; i++)
	for (int j = 1; j <= L; j++)
	if ((i + j) % 2 == 0)
		V_new[i][j] = (1 - omega) * V[i][j] + omega / 4 *
		(V[i - 1][j] + V[i + 1][j] +
		V[i][j - 1] + V[i][j + 1] +
		h * h * rho[i][j]);

	// update odd sites using updated even sites
	for (int i = 1; i <= L; i++)
	for (int j = 1; j <= L; j++)
	if ((i + j) % 2 != 0)
		V_new[i][j] = (1 - omega) * V[i][j] + omega / 4 *
		(V_new[i - 1][j] + V_new[i + 1][j] +
		V_new[i][j - 1] + V_new[i][j + 1] +
		h * h * rho[i][j]);
}

void iterate(void(*method)())
{
	clock_t t0 = clock();

	while (true) {
		method();
		++steps;
		double error = relative_error();
		if (error < accuracy)
			break;
		swap(V, V_new);         // use <algorithm> std::swap
	}
	cout << " Number of steps = " << steps << endl;

	clock_t t1 = clock();
	cout << " CPU time = " << double(t1 - t0) / CLOCKS_PER_SEC
		<< " sec" << endl;
}

int main() {

	cout << " Iterative solution of Poisson's equation\n"
		<< " ----------------------------------------\n";
	cout << " Enter number of interior points in x or y: ";
	cin >> L;

	initialize();

	cout << endl << "Rho field" << endl;
	displayField(L + 2, L + 2, rho);

	cout << endl << "V field" << endl;
	displayField(L + 2, L + 2, V);

	cout << endl << "V_new field" << endl;
	displayField(L + 2, L + 2, V_new);

	cout << " Enter desired accuracy in solution: ";
	cin >> accuracy;
	cout << " Enter 1 for Jacobi, 2 for Gauss Seidel, 3 for SOR: ";
	int choice;
	cin >> choice;

	switch (choice) {
	case 1:
		iterate(Jacobi);
		break;
	case 2:
		iterate(Gauss_Seidel);
		break;
	case 3:
		cout << " Enter overrelaxation parameter omega: ";
		cin >> omega;
		iterate(successive_over_relaxation);
		break;
	default:
		cout << " Jacobi: " << endl;
		iterate(Jacobi);
		cout << " Gauss-Seidel: " << endl;
		initialize();
		iterate(Gauss_Seidel);
		omega = 2 / (1 + 4 * atan(1.0) / double(L));
		cout << " Successive Over Relaxation with theoretical optimum omega = "
			<< omega << endl;
		initialize();
		iterate(successive_over_relaxation);
		break;
	}



	// write potential to file
	//cout << " Potential in file poisson.data" << endl;
	//ofstream date_file("poisson.data");
	//for (int i = 0; i < L + 2; i++) {
	//	double x = i * h;
	//	for (int j = 0; j < L + 2; j++) {
	//		double y = j * h;
	//		date_file << x << '\t' << y << '\t' << V[i][j] << '\n';
	//	}
	//	date_file << '\n';
	//}
	//date_file.close();

	cout << endl << "Solution:" << endl;
	displayField(L + 2, L + 2, V);


	system("Pause");
	return 0;
}