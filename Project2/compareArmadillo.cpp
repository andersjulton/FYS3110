#include "armadillo"

using namespace std;
using namespace arma;

// copy symmetric matrix A to arma mat
mat copySymMatrixToArma(double **A, int n) {
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

// test if matrix A and armadillo's matrix armaA are identical
int compareArmaMatrix(double **A, mat armaA, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (fabs(A[i][j] - armaA(i, j)) > 1e-10) {
				printf("(%d, %d) -> (%3.f != %3.f)  ", i, j, A[i][j], armaA(i, j));
				return 1;
			}
		}
	}
	return 0;
}

// test if array b and armadillo's vector armab are identical
int compareArmaVector(double *b, vec armab, int n) {
	for (int i = 0; i < n; i++) {
		if (fabs(b[i] - armab(i)) > 1e-10) {
			printf("(%d) -> (%.8f != %.8f)  \n", i, b[i], armab(i));
		return 1;
		}
	}
	return 0;
}
