#include "jacobi.h"
#include "utils.h"
#include "testLinalg.h"
#include "bisect.h"
#include "linalgUtils.h"
#include "armadillo"
#include <ctime>
#include <string>
#include <stdlib.h>

using namespace std;
using namespace arma;

void testBisect(double a, double d, int n);
void testJacobi(double a, double d, int n);
void printResult(int success, string method);

int main() {
	double a = 17;
	double d = 9;
	int n = 200;
	double **A;
	printf("Testing copySymMatrixToArma\n");
	A = createTriDiaMatrix(a, d, n);
	mat armaA = copySymMatrixToArma(A, n);
	int success = compareArmaMatrix(A, armaA, n);
	printResult(success, "copySymMatrixToArma");
	printf("Testing Jacobi\n");
	testJacobi(1, 2, 100);
	testJacobi(-4, 18, 20);
	testJacobi(8, 1, 10);
	printf("Testing Bisect\n");
	printf("n = %d\n", n);
	testBisect(a, d, n);

	double **R = createDiaMatrix(1, n);
	jacobi(A, R, n);
	double *eig = diagToVector(A, n);
	double **E = transpose(R, n);
	sortEig(eig, E, n);
	return 0;
}

void printResult(int success, string method) {
	if (success == 0) {
		printf("\t- Method %s passed the test\n\n", method.c_str());
	} else {
		printf("\t- Method %s failed the test\n\n", method.c_str());
	}

}


void testBisect(double a, double d, int n) {
	double *off = createVector(a, n);
	double *dia = createVector(d, n);
	double *eig = Gershgorin(off, dia, n);
	printf("Test eigenvalues Bisect\n");
	double **A;
	A = createTriDiaMatrix(a, d, n);
	mat armaA = copySymMatrixToArma(A, n);
	vec eigval = eig_sym(armaA);
	compareArmaVector(eig, eigval, n);

	delete[] off;
	delete[] dia;
	delete[] eig;
}

void testJacobi(double a, double d, int n) {
	double **A;
	int success;
	A = createTriDiaMatrix(a, d, n);
	mat armaA = copySymMatrixToArma(A, n);
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
