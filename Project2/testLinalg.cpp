#include "armadillo"

using namespace std;
using namespace arma;

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

// Test if Ax = lx
int testEig(double *eigval, double **eigvec, double **A, int n) {
	double left, right;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			right = eigvec[i][j]*eigval[i];
			left = 0.0;
			for (int k = 0; k < n; k++) {
				left = left + A[j][k]*eigvec[i][k];
			}
			if (fabs(right - left) > 1e-7) {
				printf("eig(%d)v(%d), %.11f != %.11f \n", i, j, right, left);
				return 1;
			}
		}
	}
	return 0;
}

// test if two arrays are equal
int testVector(double *u, double *v, int n) {
	double error;
	for (int i = 0; i < n; i++) {
		error = u[i] - v[i];
		if (fabs(error) > 1e-10) {
			printf("(%d), error = %.30f\n", i, error);
			return 1;
		}
	}
	return 0;
}


// Test orthogonality provided that A is a list of columns
int testOrthogonality(double **A, int n){
	double delta;
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			delta = 0.0;
			for (int k = 0; k < n; k++) {
				delta += A[k][i]*A[k][j];
			}
			if (j == i) {
				if (fabs(delta) < 1e-10) {
					printf("delta[%d][%d] = %.13f\n", i, j, delta);
					return 1;
				}
			} else {
				if (fabs(delta) > 1e-10){
					printf("delta[%d][%d] = %.13f\n ", i, j, delta);
					return 1;
				}
			}
		}
	}
	return 0;
}
