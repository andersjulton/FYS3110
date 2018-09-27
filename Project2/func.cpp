#include <cmath>
#include <stdlib.h>
#include "utils.h"
#include "func.h"
#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>

using namespace std;


int jacobi(double **A, double **R, int n) {
	// Setting up the eigenvector matrix
	int k, l;
	double epsilon = 1.0e-12;
	double max_number_iterations = (double) n*(double) n*(double) n;
	int iterations = 0;
  int *indexOfMax = getMaxInRow(A, n);
  k = 0;
  l = indexOfMax[0];
	double max_offdiag = maxOffDiag(A, indexOfMax, &k, &l, n);
	while ( fabs(max_offdiag) > epsilon && (double) iterations < max_number_iterations ) {
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
  if ( A[k][l] != 0.0 ) {
    double t, theta;
    theta = (A[l][l] - A[k][k])/(2*A[k][l]);
    if ( theta > 0 ) {
      t = -theta + sqrt(1.0 + theta*theta);
    } else {
      t = -theta - sqrt(1.0 + theta*theta);
    }
    c = 1/sqrt(1 + t*t);
    s = c*t;
  } else {
    c = 1.0;
    s = 0.0;
  }
  double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
  a_kk = A[k][k];
  a_ll = A[l][l];
  // changing the matrix elements with indices k and l
  A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll;
  A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll;
  A[k][l] = 0.0; // hard-coding of the zeros
  A[l][k] = 0.0;
  // and then we change the remaining elements
  for ( int i = 0; i < n; i++ ) {
    if ( i != k && i != l ) {
      a_ik = A[i][k];
      a_il = A[i][l];
      A[i][k] = c*a_ik - s*a_il;
      A[k][i] = A[i][k];
      A[i][l] = c*a_il + s*a_ik;
      A[l][i] = A[i][l];
    }
    // Finally, we compute the new eigenvectors
    r_ik = R[i][k];
    r_il = R[i][l];
    R[i][k] = c*r_ik - s*r_il;
    R[i][l] = c*r_il + s*r_ik;
  }
  return;
}

double maxOffDiag(double **A, int *indexOfMax, int *k, int *l, int n) {
  for (int i = 0; i < n; i++) {
    if (fabs(A[i][indexOfMax[i]]) > fabs(A[*k][*l])) {
      *k = i;
      *l = indexOfMax[i]; 
    }
  }
  return A[*k][*l];
}

void updateMaxInRow(double **A, int *indexOfMax, int k, int l, int n) {
  // update index of largest element in row l & k
  for (int i = 0; i < n; i++) {
    if (fabs(A[l][i]) > fabs(A[l][indexOfMax[l]]) && i != l) {
      indexOfMax[l] = i; 
    } 
    if (fabs(A[k][i]) > fabs(A[k][indexOfMax[k]]) && i != k) {
      indexOfMax[k] = i; 
    } // see if update in column made changes in max in row i
    if (indexOfMax[i] == k || indexOfMax[i] == l) {
      // Update row if l or k justed to be max in row 
      if (i != k && i != l) {
        // row l and k are regardless updated
        for (int j = 0; j < n; j++) {
          if (fabs(A[i][j]) > fabs(A[i][indexOfMax[i]]) && i != j) {
            indexOfMax[i] = j;
          }
        }
      }
    } else {
      // upate max in row if changed col-value is larger
      if ((fabs(A[i][l]) > fabs(A[i][indexOfMax[i]])) && i != l) {
        indexOfMax[i] = l;
      }
      if ((fabs(A[i][k]) > fabs(A[i][indexOfMax[i]])) && i != k) {
        indexOfMax[i] = k;
      }
    }
  }
}

// get position of the largest value in row when i != j
int* getMaxInRow(double **A, int n) {
  int *indexOfMax = new int[n];
  for (int i = 0; i < n; i++) {
    if (i > 0) {
      indexOfMax[i] = 0;
    } else {
      indexOfMax[0] = 1;
    }
    for (int j = 0; j < n; j++) {
      if (fabs(A[i][j]) > fabs(A[i][indexOfMax[i]]) && i != j) {
        indexOfMax[i] = j;
      }
    }
  }
  return indexOfMax;
}
