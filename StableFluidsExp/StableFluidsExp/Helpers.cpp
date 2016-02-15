#include "Helpers.h"

void displayField(int row, int col, float* field){
	cout << std::setprecision(1) << setiosflags(ios::fixed);
	for (int i = 0; i != row; i++){
		for (int j = 0; j != col; j++){
			cout << std::setw(7) << field[i * col + j] << ", ";
			if (j == col - 1)
				cout << endl;
		}
	}
}

void displayVectorField(int row, int col, float* u, float* v){
	cout << std::setprecision(2) << setiosflags(ios::fixed);
	for (int i = 0; i != row; i++){
		for (int j = 0; j != col; j++){
			cout << std::setw(3) << "v(" << u[i * col + j] << "," << v[i * col + j] << ") ";
			if (j == col - 1)
				cout << endl;
		}
	}
}

void displayFieldInv(int row, int col, float * field){
	cout << std::setprecision(1) << setiosflags(ios::fixed);
	for (int i = row - 1; i >= 0; i--){
		for (int j = 0; j != col; j++){
			cout << std::setw(7) << field[i * col + j] << ", ";
			if (j == col - 1)
				cout << endl;
		}
	}
}

void displayVectorFieldInv(int row, int col, float * u, float * v){
	cout << std::setprecision(2) << setiosflags(ios::fixed);
	for (int i = row - 1; i >= 0; i--){
		for (int j = 0; j != col; j++){
			cout << std::setw(3) << "v(" << u[i * col + j] << "," << v[i * col + j] << ") ";
			if (j == col - 1)
				cout << endl;
		}
	}
}

void computeCurls_uniform(int dim, float * w, float * u, float * v)
{
	float h = 1.0f / dim;

	int N = dim - 3;

	for (int i = N + 1; i >= 1; i--)
	{
		for (int j = 1; j <= N + 1; j++)
		{
			float du, dv, curl;
			float uTop, vLeft, uBottom, vRight;
			uBottom = (u[IX(i - 1, j - 1)] + u[IX(i - 1, j)]) * 0.5f;
			uTop = (u[IX(i, j)] + u[IX(i, j - 1)]) * 0.5f;
			vLeft = (v[IX(i, j - 1)] + v[IX(i - 1, j - 1)]) * 0.5f;
			vRight = (v[IX(i - 1, j)] + v[IX(i, j)]) * 0.5f;
			du = (uBottom - uTop) / h;
			dv = (vRight - vLeft) / h;
			curl = dv - du;
			w[IXCURL(i, j)] = curl;
		}
	}
}