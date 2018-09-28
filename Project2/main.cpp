#include "jacobi.h"
#include "utils.h"
#include <string>

using namespace std;

int main() {
	int n = 200;
	int iterations;
	double rho_0, rho_N, rho_i, h;
	double **A, **R;
	double *a;
	a = createVector(0, n);
	R = createDiaMatrix(1, n);
	rho_N = 10.0;
	rho_0 = 0;
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
		cout << a[j] << '\n';
	}
	deleteMatrix(A, n);
	deleteMatrix(R, n);
	delete [] a;
	return 0;
}
