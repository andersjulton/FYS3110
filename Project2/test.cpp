#include "jacobi_bisect.h"
#include "utils.h"
#include "compareArmadillo.h"
#include "linalgUtils.h"
#include "armadillo"
#include <string>
#include <stdlib.h>
#include <cmath>

using namespace std;
using namespace arma;

void testBisect(double a, double d, int n);
void testJacobi(double a, double d, int n);
void printResult(int success, string method);
void iterToFileJacobi(double a, double d);
void testRotate(double a, double d);
void testFindMax(int n);
void testCreateTriDiaMatrix();

int main() {
	printf("\n------------------------------------------------------------------------\n");
	testCreateTriDiaMatrix();

	double a = -1;
	double d = 2;
	int n = 20;

	printf("\n------------------------------------------------------------------------\n");
	testFindMax(n);

	printf("\n------------------------------------------------------------------------\n");
	testRotate(a, d);


	printf("\n------------------------------------------------------------------------\n");
	testJacobi(1, 2, 100);

	printf("\n------------------------------------------------------------------------\n");
	testBisect(a, d, n);

	//iterToFileJacobi(-1, 2);
	return 0;
}

void testCreateTriDiaMatrix() {
	int n = 20;
	double pi = 3.14159265359;
	double h = 1.0/n;
	double a =  -1.0/(h*h);;
	double d = 2.0/(h*h);
	printf("Testing transfer of values from created matrix A to a armadillo matrix\n");
	double **A = createTriDiaMatrix(a, d, n);
	mat armaA = copySymMatrixToArma(A, n);
	int success = compareArmaMatrix(A, armaA, n);
	printResult(success, "copySymMatrixToArma");

	printf("Testing if created matrix gives expected eigenvalues\n");
	vec eigval = eig_sym(armaA);
	double *expected = new double[n];
	for (int i = 0; i < n; i++) {
		expected[i] = d + 2*a*cos((i+1)*pi/(n + 1));
	}
	success = compareArmaVector(expected, eigval, n);

	printResult(success, "createTriDiaMatrix");
	deleteMatrix(A, n);
	delete[] expected;
}

void printResult(int success, string method) {
	if (success == 0) {
		printf("\t- Method %s passed the test\n", method.c_str());
	} else {
		printf("\t- Method %s failed the test\n", method.c_str());
	}
}

void testBisect(double a, double d, int n) {
	int success;
	double *off = createVector(a, n);
	double *dia = createVector(d, n);
	double *eig = bisect(off, dia, n, 1e-13);
	double **A = createTriDiaMatrix(a, d, n);
	mat armaA = copySymMatrixToArma(A, n);
	vec eigval = eig_sym(armaA);

	printf("Checking if Bisect returns expected eigenvalues\n");
	success = compareArmaVector(eig, eigval, n);
	printResult(success, "bisect");

	deleteMatrix(A, n);
	delete[] off;
	delete[] dia;
	delete[] eig;
}

void testJacobi(double a, double d, int n) {
	int success;
	double **A = createTriDiaMatrix(a, d, n);
	mat armaA = copySymMatrixToArma(A, n);
	double **R = createDiaMatrix(1, n);

	jacobi(A, R, n, 1e-13);
	double *eig = diagToVector(A, n);
	double **E = transpose(R, n);
	sortEig(eig, E, n);
	vec eigval = eig_sym(armaA);

	printf("Checking if the Jacobi algorithm returns expected eigenvalues\n");
	success = compareArmaVector(eig, eigval, n);
	printResult(success, "Jacobi");

	printf("Testing Orthogonality of eigenvectors returned by Jacobi\n");
	success = testOrthogonality(E, n);
	printResult(success, "Jacobi");

	printf("Testing Ax = lx for the returned eigenvalues and eigenvectors\n");
	deleteMatrix(A, n);		// A has been changed
	A = createTriDiaMatrix(a, d, n);
	success = testEig(eig, E, A, n);
	printResult(success, "Jacobi");

	delete[] eig;
	deleteMatrix(E, n);
	deleteMatrix(A, n);
	deleteMatrix(R, n);
}

void testFindMax(int n) {
	double **A = createTriDiaMatrix(1.0, 60.0, n);
	int *maxPos = getMaxInRow(A, n);
	int k = 0;
	int l = maxPos[k];
	printf("Testing if Jacobi's supporting functions distinguishes between diagonal and off-diagonal elements\n");
	double max = maxOffDiag(A, maxPos, &k, &l, n);
	int success = 0;
	if (fabs(60.0 - max) < 1e-12) {
		success = 1;
	}
	printResult(success, "getMaxInRow");

	printf("Testing if Jacobi's supporting functions finds largest off-diagonal element\n");
	delete[] maxPos;
	int i = n/3;
	int j = i + 2;
	double myMax = 3.14;
	A[i][j] = myMax;
	maxPos = getMaxInRow(A, n);
	max = maxOffDiag(A, maxPos, &k, &l, n);
	success = 0;
	if (fabs(myMax - max) > 1e-12) {
		success = 1;
	}
	printResult(success, "getMaxInRow");

	printf("Testing if Jacobi's supporting functions finds the correct coordinates of value max\n");
	success = (i - k) + (j - l);
	printResult(success, "maxOffDiag");

	printf("Testing if Jacobi's supporting functions updates max-coordinates correctly\n");
	myMax = 42.0;
	i = k;
	j = i + 1;
	A[i][j] = myMax;
	updateMaxInRow(A, maxPos, k, l, n);
	max = maxOffDiag(A, maxPos, &k, &l, n);
	success = 0;
	if (fabs(myMax - max) > 1e-12) {
		success = 1;
	}
	printResult(success, "updateMaxInRow (in row)");
	myMax = 9000.0;
	i = k - 1;
	j = l;
	A[i][j] = myMax;
	updateMaxInRow(A, maxPos, k, l, n);
	max = maxOffDiag(A, maxPos, &k, &l, n);
	success = 0;
	if (fabs(myMax - max) > 1e-12) {
		success = 1;
	}
	printResult(success, "updateMaxInRow (in column)");
	deleteMatrix(A, n);
	delete[] maxPos;
}

void testRotate(double a, double d) {
	int n = 10;
	double **A = createTriDiaMatrix(a, d, n);
	double **R = createDiaMatrix(1, n);
	double **E;
	double maxOffdiag;
	int k, l;
	int *indexOfMax = getMaxInRow(A, n);
	k = 0;
	l = indexOfMax[0];
	int success = 0;
	for (int i = 0; i < 5; i++) {
		maxOffdiag = maxOffDiag(A, indexOfMax, &k, &l, n);
		rotate(A, R, k, l, n);
		updateMaxInRow(A, indexOfMax, k, l, n);

		// test
		E = transpose(R, n);
		success += testOrthogonality(E, n);
		deleteMatrix(E, n);
	}
	printf("Testing if Orthogonality is preseved after every five rotation in the Jacobi algorithm\n"); 
	printResult(success, "rotate");
	
	delete[] indexOfMax;
	deleteMatrix(A, n);
	deleteMatrix(R, n);
}


void iterToFileJacobi(double a, double d) {
	int *n_vector = new int[40];
	int *iter_vector = new int[40];
	int n = 2;
	double **R, **A;
	for (int i = 0; i < 10; i++) {
		n_vector[i] = n;
		A = createTriDiaMatrix(a, d, n);
		R = createDiaMatrix(1, n);
		iter_vector[i] = jacobi(A, R, n, 1e-12);
		deleteMatrix(A, n);
		deleteMatrix(R, n);
		n += 2;
	}
	n = 40;
	for (int i = 0; i < 30; i++) {
		n_vector[i+10] =  n;
		A = createTriDiaMatrix(a, d, n);
		R = createDiaMatrix(1, n);
		iter_vector[i+10] = jacobi(A, R, n, 1e-12);
		deleteMatrix(A, n);
		deleteMatrix(R, n);
		n += 20;
	}
	intArrayToFile(n_vector, 40, "n");
	intArrayToFile(iter_vector, 40, "iter");
	delete[] n_vector;
	delete[] iter_vector;
}
