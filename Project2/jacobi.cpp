#include <cmath>
#include <algorithm>
#include "utils.h"
#include "jacobi.h"


int jacobi(double **A, double **R, int n, double epsilon) {
	int k, l;
	int max_number_iterations = n*n*n;
	int iterations = 0;
	int *indexOfMax = getMaxInRow(A, n);
	k = 0;
	l = indexOfMax[0];
	double max_offdiag = maxOffDiag(A, indexOfMax, &k, &l, n);
	while (fabs(max_offdiag) > epsilon && iterations < max_number_iterations ) {
		rotate(A, R, k, l, n);
		updateMaxInRow(A, indexOfMax, k, l, n);
		max_offdiag = maxOffDiag(A, indexOfMax, &k, &l, n);
		iterations++;
	}
	delete[] indexOfMax;
	return iterations;
}

void rotate(double **A, double **R, int k, int l, int n) {
	double s, c;
	// Calculate phi
	if (A[k][l] != 0.0) {
		double t, phi;
		phi = (A[l][l] - A[k][k])/(2*A[k][l]);
		if ( phi > 0 ) {
			t = -phi + sqrt(1.0 + phi*phi);
		} else {
			t = -phi - sqrt(1.0 + phi*phi);
		}
		c = 1/sqrt(1 + t*t);
		s = c*t;
	} else {             // This should NEVER happen
		c = 1.0;
		s = 0.0;
	}
	double A_kk, A_ll, A_ik, A_il, R_ik, R_il;
	A_kk = A[k][k];
	A_ll = A[l][l];
	// Rotate matrix and update matrix elements
	A[k][k] = c*c*A_kk - 2.0*c*s*A[k][l] + s*s*A_ll;
	A[l][l] = s*s*A_kk + 2.0*c*s*A[k][l] + c*c*A_ll;
	A[k][l] = 0.0;
	A[l][k] = 0.0;
	for ( int i = 0; i < n; i++ ) {
		if ( i != k && i != l ) {
			A_ik = A[i][k];
			A_il = A[i][l];
			A[i][k] = c*A_ik - s*A_il;
			A[k][i] = A[i][k];
			A[i][l] = c*A_il + s*A_ik;
			A[l][i] = A[i][l];
		}
		// Computing new eigenvectors
		R_ik = R[i][k];
		R_il = R[i][l];
		R[i][k] = c*R_ik - s*R_il;
		R[i][l] = c*R_il + s*R_ik;
	}
}

double maxOffDiag(double **A, int *indexOfMax, int *k, int *l, int n) {
	for (int i = 0; i < n-1; i++) {
		if (fabs(A[i][indexOfMax[i]]) > fabs(A[*k][*l])) {
			*k = i;
			*l = indexOfMax[i];
		}
	}
	return A[*k][*l];
}

void updateMaxInRow(double **A, int *indexOfMax, int k, int l, int n) {
	// see if max in column was k or l
	for (int i = 0; i < (n-1); i++) {
		if (indexOfMax[i] == k || indexOfMax[i] == l) {
			// Update row if l or k used to be max in row
			for (int j = i+1; j < n; j++) {
				if (fabs(A[i][j]) > fabs(A[i][indexOfMax[i]])) {
					indexOfMax[i] = j;
				}
			}
		}
	}
	// see if new element in column k or l > max
	for (int i = 0; i < k; i++) {
		if (std::max(fabs(A[i][k]), fabs(A[i][k])) > fabs(A[i][indexOfMax[i]])) {
			if ( (fabs(A[i][k]) > fabs(A[i][l])) ) {
				indexOfMax[i] = k;
			} else {
				indexOfMax[i] = l;
			}
		} //else {printf("This happens\n");}
	}
	//printf("%d\n",(l-k) );
	// k < l
	for (int i = k; i < l; i++) {
		if ((fabs(A[i][l]) > fabs(A[i][indexOfMax[i]]))) {
			indexOfMax[i] = l;
		}
	}
}

// get position of the largest value in row of the upper triagonal part when i != j
int* getMaxInRow(double **A, int n) {
	int *indexOfMax = new int[n-1];
  	for (int i = 0; i < (n - 1); i++) {
  		indexOfMax[i] = i + 1;
    	for (int j = i + 2; j < n; j++) {
    		if (fabs(A[i][j]) > fabs(A[i][indexOfMax[i]])) {
				indexOfMax[i] = j;
			}
    	}
  	}
	return indexOfMax;
}
