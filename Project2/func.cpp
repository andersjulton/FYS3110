#include <cmath>
#include <stdlib.h>
#include "utils.h"
#include "func.h"
#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>

using namespace std;

int main () {
	double **A;
	int n = 5;
	int iterations;
	A = createTriDiaMatrix(1, 2, n);
	iterations = jacobi(A, n);
	std::cout << A[0][0] << '\n';
	for (int i = 0; i < n; i++) {
		std::cout << A[i][i] << '\n';
	}
	return 0;
}
int jacobi(double **A, int n) {
	// Setting up the eigenvector matrix
	int k, l;
	double **R;
	R = createDiaMatrix(n, 1);
	double epsilon = 1.0e-8;
	double max_number_iterations = (double) n*(double) n*(double) n;
	int iterations = 0;
	double max_offdiag = maxOffDiag(A, &k, &l, n);
	while ( fabs(max_offdiag) > epsilon && (double) iterations < max_number_iterations ) {
		max_offdiag = maxOffDiag(A, &k, &l, n);
		rotate(A, R, k, l, n);
		iterations++;
	}
	return iterations;
}

void rotate(double **A, double **R, int k, int l, int n) {
  double s, c;
  if ( A[k][l] != 0.0 ) {
    double t, theta;
    theta = (A[l][l] - A[k][k])/(2*A[k][l]);
    if ( theta > 0 ) {
      t = 1.0/(theta + sqrt(1.0 + theta*theta));
    }
    else {
      t = -1.0/(theta + sqrt(1.0 + theta*theta));
    }
    c = 1/sqrt(1 + t*t);
    s = c*t;
  }
  else {
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

double maxOffDiag(double **A, int *k, int *l, int n) {
  double max = 0.0;
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      if (fabs(A[i][j]) > max) {
        max = fabs(A[i][j]);
        *l = i;
        *k = j;
      }
    }
  }
  return max;
}

void bisect(double *a, double *b, int n) { // should it return the eigenvalues?
	// Get interval by the Gershgorin circle theorem
	double *interval;
	int xmin, xmax;
	interval = getInterval(a, b, n);
	xmin = interval[0];
	xmax = interval[1];
	// Calculate number of eigenvalues
	int smin, smax, nlambda;
	smin = getSignChange(a, b, n, xmin);
	smax = getSignChange(a, b, n, xmax);
	nlambda = smin - smax;
	// Calculate eigenvalues
	double eps = 10e-8;
	double q0, q1, x1, xu, x0;
	x0 = xmax;
	int z = 0;		// what is z? A counter?
	int m1 = 0;
	int u;
	double *x, *wu;
	x = createVector(xmax, nlambda);
	wu = createVector(xmin, nlambda);
	for (int k = nlambda - 1; k > -1; k--) {
		xu = xmin;
		for (int i = k; i > -1; i--) {
			if (xu < wu[i]) {
				xu = wu[i];
			}
			if (x0 > x[k]) {
				x0 = x[k];
			}
			while ((x0 - xu) > eps) {
				z += 1;			// counter?
				u = 0;
				q0 = 1;			// variable init. too many times. Should it change?
				for (int j = 0; j < n; j++) {
					if (q0 == 0) {	// why is this in the loop? q0 is always 1. This is NEVER true
						q0 = eps;
					}
					q1 = a[i] - x1 - b[i]*b[i]/q0;  // q0 is always one
					if (q1 < 0) {
						u += 1;
					}
					if (u < k) {
						if (u < m1) {	// so all times but 1?
							xu = wu[m1];
							x1 = xu;
						} else {
							xu = wu[u + 1];
							x1 = xu;
							if (x[u] > x1) {
								x[u] = x1;
							}
						}
					} else {
						x0 = x1;
					}
				}
			}
		}
		x[k] = (x0 + xu)/2.0;		// should x be returned?
	}
	delete[] x; 					// should x be returned?
	delete[] wu;
}

// Find max value of interval by Gershgorin’s Theorem
double *getInterval(double *a, double *b, int n) {
  double xmax, xmin, h;
  double *interval;
  interval = createVector(0, 2);
  xmin = a[n-1] - abs(b[n-1]);
  xmax = a[n-1] + abs(b[n-1]);
  for (int i = (n - 2); i > -1; i--) {
      h = abs(b[i]) + abs(b[i+1]);
      if ((a[i] + h) > xmax) {
        xmax = a[i] + h;
      }
      if ((a[i] - h) < xmin) {
        xmin = a[i] - h;
      }
  }
  interval[0] = xmin;
  interval[1] = xmax;
  return interval;
}

int getSignChange(double *a, double *b, int n, double lambda) {
  double *p, *s;
  int sum = 0;			// sum or count?
  s = createVector(0, n+1);
  p = createVector(0, n+1);
  p[0] = 1; p[1] = a[0] - lambda;
  s[0] = 1;
  if (p[1] < 0) {
    s[1] = -1;
  } else if (p[1] > 0 ) {
    s[1] = 1;
  } else {
    s[1] = -1*s[0];
  }
  for (int i = 2; i < n+1; i++) {
      p[i] = (a[0] - lambda)*p[i-1] - p[i-2];
      if (p[i] < 0) {
        s[i] = -1;
      } else if (p[i] > 0) {
        s[i] = 1;
      } else {
        s[i] = -1*s[i-1];
      }
  }
  for (int i = 0; i < n; i++) {
    if(((s[i] > 0) && (s[i+1] > 0)) || ((s[i] < 0) && (s[i+1] < 0))) {
      sum += 1;
    }
  }
  delete[] s;
  delete[] p;
  return sum;
}
