#include <cmath>
#include <stdlib.h>
#include <algorithm>
#include <iostream>

double **createMatrix(int m, int n);
void deleteMatrix(double **mat, int n);
double **createDiaMatrix(double value, int n);
double **createTriDiaMatrix(double off_value, double d_value, int n);
double *createVector(double value, int n);
int jacobi(double **A, double **R, int n);
void rotate(double **A, double **R, int k, int l, int n);
double maxOffDiag(double **A, int *k, int *l, int n);
void bisect(double *a, double *b, int n);
int getSignChange(double *a, double *b, int n, double lambda);
double *getInterval(double *a, double *b, int n);
void *vectorMatrixMultiplication(double *a, double **A, double *u, int n);
double *createEVector(int index, int n);
double innerProduct(double *a, double *b, int n);
void lanczos(double **A, int n);
void scalarVector(double *v, double *u, double value, int n);



int main() {
    int n = 5;
    double **A;
    A = createTriDiaMatrix(-1, 4, n);
    lanczos(A, n);
  /*int n = 200;
  int iterations;
  double rho_0, rho_N, rho_i, h;
  double **A, **R;
  double *a;
  a = createVector(0, n);
  R = createDiaMatrix(1, n);
  rho_0 = h; rho_N = 10.0;
  h = (rho_N - rho_0)/n;

  A = createMatrix(n, n);
  A[0][0] = 2.0/(h*h) + rho_0*rho_0;
  A[0][1] = -1.0/(h*h);
  for (int i = 1; i < (n - 1); i++) {
    rho_i = (i+1)*h;
    A[i][i-1] = -1.0/(h*h);
    A[i][i] = 2.0/(h*h) + rho_i*rho_i;
    A[i][i+1] = -1.0/(h*h);
  }
  rho_i = (n - 1);
  A[n-1][n-2] = -1.0/(h*h);
  A[n-1][n-1] = 2.0/(h*h) + rho_i*rho_i;

  iterations = jacobi(A, R, n);
  for (int j = 0; j < n; j++) {
    a[j] = A[j][j];
  }
  std::sort(a, a+n);
  for (int j = 0; j < 10; j++) {
    std::cout << a[j] << '\n';
  }
  deleteMatrix(A, n);
  deleteMatrix(R, n);
  delete [] a;*/
  return 0;
}

//  Allocating space for a m x n matrix and fill elements with 0
double **createMatrix(int m, int n) {
	double **mat;
	mat = new double*[m];
	for(int i = 0; i < m; i++) {
		mat[i] = new double[n];
		for(int j = 0; j < n; j++) {
			mat[i][j] = 0.0;
		}
	}
	return mat;
}

// delete 2D-array
void deleteMatrix(double **mat, int n) {
	for (int i = 0; i < n; i++) {
		delete[] mat[i];
	}
	delete[] mat;
}

double **createDiaMatrix(double value, int n) {
  double **mat;
  mat = new double*[n];
  for(int i = 0; i < n; i++) {
    mat[i] = new double[n];
    for(int j = 0; j < n; j++) {
      if ( i == j ) {
        mat[i][j] = value;
      }
      else {
        mat[i][j] = 0.0;
      }
    }
  }
  return mat;
}

// Create a tridiagonal matrix
double **createTriDiaMatrix(double off_value, double d_value, int n) {
	// two seperate loops to avoid if-tests
	double **mat;
	mat = createMatrix(n, n);
	mat[0][0] = d_value;
	mat[0][1] = off_value;
	for (int i = 1; i < (n - 1); i++) {
		mat[i][i-1] = off_value;
		mat[i][i] = d_value;
		mat[i][i+1] = off_value;
	}
	mat[n-1][n-2] = off_value;
	mat[n-1][n-1] = d_value;
	return mat;
}

double *createVector(double value, int n) {
	//  Allocating space for a vector and fill elements with a constant value
	double *vector;
	vector = new double[n];
	for (int i = 0; i < n; i++) {
		vector[i] = value;
	}
	return vector;
}

int jacobi(double **A, double **R, int n) {
	// Setting up the eigenvector matrix
	int k, l;
	double epsilon = 1.0e-12;
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
      t = -theta + sqrt(1.0 + theta*theta);
    }
    else {
      t = -theta - sqrt(1.0 + theta*theta);
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

void bisect(double *a, double *b, int n) {
  // Get interval by the Gershgorin circle theorem
  double *interval;
  int xmin, xmax;
  interval = getInterval(a, b, n);
  xmin = interval[0]; xmax = interval[1];
  // Calculate number of eigenvalues
  int smin, smax, nlambda;
  smin = getSignChange(a, b, n, xmin);
  smax = getSignChange(a, b, n, xmax);
  nlambda = smin - smax;
  // Calculate eigenvalues
  double eps = 10e-12;
  double q0, q1, x1, xu, x0;
  x0 = xmax;
  int z = 0;
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
        z += 1;
        u = 0; q0 = 1;
        for (int j = 0; j < n; j++) {
          if (q0 == 0) {
            q0 = eps;
          }
          q1 = a[i] - x1 - b[i]*b[i]/q0;
          if (q1 < 0) {
            u += 1;
          }
          if (u < k) {
             if (u < m1) {
              xu = wu[m1];
              x1 = xu;
            }
            else {
              xu = wu[u + 1];
              x1 = xu;
              if (x[u] > x1) {
                x[u] = x1;
              }
            }
          }
          else {
            x0 = x1;
          }
        }
      }
    }
    x[k] = (x0 + xu)/2.0;
  }
}

double *getInterval(double *a, double *b, int n) {
  double xmax, xmin, h;
  double *interval;
  interval = createVector(0, 2);
  xmin = a[n-1] - abs(b[n-1]);
  xmax = a[n-1] + abs(b[n-1]);
  // Find max value of interval by Gershgorin’s Theorem
  for (int i = n - 2; i > -1; i--) {
      h = abs(b[i]) + abs(b[i+1]);
      if ((a[i] + h) > xmax) {
        xmax = a[i] + h;
      }
      if ((a[i] - h) < xmin) {
        xmin = a[i] - h;
      }
  }
  interval[0] = xmin; interval[1] = xmax;
  return interval;
}

int getSignChange(double *a, double *b, int n, double lambda) {
  double *p, *s;
  int sum = 0;
  s = createVector(0, n+1);
  p = createVector(0, n+1);
  p[0] = 1; p[1] = a[0] - lambda;
  s[0] = 1;
  if (p[1] < 0) {
    s[1] = -1;
  }
  else if (p[1] > 0 ) {
    s[1] = 1;
  }
  else if (p[1] == 0) {
    s[1] = -1*s[0];
  }
  for (int i = 2; i < n+1; i++) {
      p[i] = (a[0] - lambda)*p[i-1] - p[i-2];
      if (p[i] < 0) {
        s[i] = -1;
      }
      else if (p[i] > 0 ) {
        s[i] = 1;
      }
      else if (p[i] == 0) {
        s[i] = -1*s[i-1];
      }
  }
  for (int j = 0; j < n; j++) {
    if(((s[j] > 0) && (s[j+1] > 0)) || ((s[j] < 0) && (s[j+1] < 0))) {
      sum += 1;
    }
  }
  return sum;
}

void lanczos(double **A, int n) {
    double *a, *v, *w_, *w, *b;
    a = createVector(0, n);
    v = createEVector(0, n);
    w_ = createVector(0, n);
    w = createVector(0, n);
    b = createVector(0, n);
    vectorMatrixMultiplication(v, A, w_, n);
    a[0] = innerProduct(w_, v, n);
    for (int i = 0; i < n; i++) {
        w[i] = w_[i] - a[0]*v[i];
    }
    for (int i = 1; i < n; i++) {
        b[i] = fabs(innerProduct(w, w, n));
        if (b[i] != 0) {
            scalarVector(v, w, 1.0/b[i], n);
        } else {
            v = createEVector(i, n);
        }
        vectorMatrixMultiplication(v, A, w_, n);
        a[i] = innerProduct(w_, v, n);
        for (int j = 0; j < n; j++) {
            w[j] = w_[j] - a[i]*v[j];
        }
    }
    for (int i = 0; i < n; i++) {
        std::cout << a[i] << '\n';
    }
}

double *createEVector(int index, int n) {
	//  Allocating space for a vector and fill elements with a constant value
	double *vector;
	vector = new double[n];
	for (int i = 0; i < n; i++) {
		vector[i] = 0;
	}
    vector[index] = 1;
	return vector;
}

void *vectorMatrixMultiplication(double *a, double **A, double *u, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            u[i] += A[i][j]*a[j];
        }
    }
}
double innerProduct(double *a, double *b, int n) {
    double dot_prod = 0.0;
    for (int i = 0; i < n; i++) {
        dot_prod += a[i]*b[i];
    }
    return dot_prod;
}
void scalarVector(double *v, double *u, double value, int n) {
    for (int i = 0; i < n; i++) {
        u[i] = value*v[i];
    }
}
