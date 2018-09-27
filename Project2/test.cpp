#include <fstream>
#include "func.h"
#include "utils.h"
#include "testFunc.h"
#include "bisect.h"
#include "linalgUtils.h"
#include "armadillo"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <string>

using namespace std;
using namespace arma;

void testBisect(double a, double d, int n);
void testJacobi(double a, double d, int n);

int main() {
	/*printf("Testing Jacobi\n");
	testJacobi(1, 2, 100);
	testJacobi(-4, 18, 20);
	testJacobi(8, 1, 10);*/
	printf("Testing Bisect\n");
	double a = -7;
	double d = 15;
	int n = 300;
	printf("n = %d\n", n);
	testBisect(a, d, n);

	double **A;
	A = createTriDiaMatrix(a, d, n);
	double **R = createDiaMatrix(1, n);
	jacobi(A, R, n);
	double *eig = diagToVector(A, n);
	double **E = transpose(R, n);
	sortEig(eig, E, n);
	/*printf("\nJacobi\n");
	for (int i = 0; i < n; i++) {
		printf("%.7f\n", eig[i]);
	}*/

	return 0;
}


void testBisect(double a, double d, int n) {
	double *off = createVector(a, n);
	double *dia = createVector(d, n);
	double *eig = Gershgorin(off, dia, n);
	printf("Test eigenvalues Bisect\n");
	double **A;
	A = createTriDiaMatrix(a, d, n);
	mat armaA = copyMatrixToArma(A, n);
	vec eigval = eig_sym(armaA);
	compareArmaVector(eig, eigval, n);


	delete[] off;
	delete[] dia;
	delete[] eig;
}

void testJacobi(double a, double d, int n) {
	double **A;
	A = createTriDiaMatrix(a, d, n);
	mat armaA = copyMatrixToArma(A, n);
	compareArmaMatrix(A, armaA, n);
	// check eigenvalues
	double **R = createDiaMatrix(1, n);
	jacobi(A, R, n);
	double *eig = diagToVector(A, n);
	double **E = transpose(R, n);
	sortEig(eig, E, n);
	vec eigval = eig_sym(armaA);

	compareArmaVector(eig, eigval, n);
	// check eigenvectors
	printf("Testing Orthogonality\n");
	testOrthogonality(E, n);

	printf("Testing AE = lamdaE\n");
	deleteMatrix(A, n);
	A = createTriDiaMatrix(a, d, n);
	testEig(eig, E, A, n);

	deleteMatrix(E, n);
	deleteMatrix(A, n);
	deleteMatrix(R, n);
	printf("\n");
}
