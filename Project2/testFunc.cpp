#include <fstream>
#include "func.h"
#include "utils.h"
#include "armadillo"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <string>

using namespace std;
using namespace arma;

// copy symmetric matrix A to arma mat
mat copyMatrixToArma(double **A, int n) {
	mat armaA(n, n, fill::zeros);

	for (int i = 0; i < n; i++) {
		armaA(i, i) = A[i][i];
		for (int j = i+1; j < n; j++) {
			armaA(i, j) = A[i][j];
			armaA(j, i) = A[i][j];
		}
	}
	return armaA;
}

int compareArmaMatrix(double **A, mat armaA, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			//printf(" (%3.f, %3.f)  ", A[i][j], armaA(i, j));
			if (fabs(A[i][j] - armaA(i, j)) > 1e-10) {
				printf("(%d, %d) -> (%3.f != %3.f)  ", i, j, A[i][j], armaA(i, j));
				return 1;
			}
		}
		//printf("\n");
	}
	return 0;
}

int compareArmaVector(double *b, vec armab, int n) {
	for (int i = 0; i < n; i++) {
		//printf(" (%3.f, %3.f)  ", b[i], armab(i));
		if (fabs(b[i] - armab(i)) > 1e-10) {
			printf("(%d) -> (%.8f != %.8f)  \n", i, b[i], armab(i));
		return 1;
		}
	}
	//printf("\n");
	return 0;
}
