/*
 * Helpers.cpp
 * Xue Cui
 * 15th Mar, 2016
 * Desc:
 *
*/
#include "Helpers.h"

void displayField(int row, int col, cpt::Matrix<double, 2> field)
{
	cout << std::setprecision(4) << setiosflags(ios::fixed);
	for (int i = 0; i != row; i++){
		for (int j = 0; j != col; j++){
			cout << std::setw(7) << field[i][j] << ", ";
			if (j == col - 1)
				cout << endl;
		}
	}
}