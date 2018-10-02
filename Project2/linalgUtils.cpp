#include <cmath>
#include <string>

using namespace std;

/* version of insertion sort for sorting a list of eigenvalues 
their corresponding eigenvectors will be placed accordingly
eigenvec must be a list of columns where each colomn is an eigenvector*/
void sortEig(double *eigval, double **eigvec, int n) {
	int j;
	double current, past;
	double *currentVec, *pastVec;
	for (int i = 1; i < n; i++) {
		j = i - 1;
		current = eigval[i];
		currentVec = eigvec[i];
		past = eigval[j];
		pastVec = eigvec[j];
		while (current < past) {
			eigval[j+1] = past;
			eigvec[j+1] = pastVec;
			j--;
			if (j < 0) {
				break;
			}
			past = eigval[j];
			pastVec = eigvec[j];
		}
		eigval[j+1] = current;
		eigvec[j+1] = currentVec;
	}
}

// transpose n x n matrix 
double** transpose(double **A, int n) {
	double **transA = new double*[n];
	for (int i = 0; i < n; i++) {
		transA[i] = new double[n];
		for (int j = 0; j < n; j++) {
			transA[i][j] = A[j][i];
		}
	}
	return transA;
}

// Normalize a vector
void normalize(double *v, double *u, int n) {
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += v[i];
	}
	for (int j = 0; j < n; j++) {
		u[j] = v[j]/sum;
	}
}

// Exctract an eigenvector from matrix
void extractEigenVec(double **A, double *u, int index, int n) {
	for (int i = 0; i < n; i++) {
		u[i] = A[index][i];
	}
}

// Represent diagonal matrix A as a vector
double* diagToVector(double **A, int n) {
	double *a = new double[n];
	for (int i = 0; i < n; i++) {
		a[i] = A[i][i];
	}
	return a;
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
