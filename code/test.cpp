#include <fstream>
#include "func.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <string>

int main() {
	int n;
	double *myArray, *otherArray;
	double eps = 1e-16;

	printf("\nTESTING\n");

	// testing *double createVector(double value, int n)
	printf("\n*double createVector(double value, int n)\n");
	n = 30;
	bool success = true;
	double value = 5.5;

	myArray = createVector(value, n);

	for (int i = 0; i < n; i++) {
		if (fabs(myArray[i] - value) > eps) {
			success = false;
			break;
		}
	}
	if (success) {
		printf("\t- PASSED : elements in array matched expected value\n\n");
	} else {
		printf("\t- FAILED : elements in array did not match expected value\n\n");
	}
	delete[] myArray;

	// testing double maxError(double *u, double *v, int n)
	printf("\ndouble maxError(double *u, double *v, int n)\n");
	myArray = createVector(1, n);
	otherArray = new double[n];
	success = true;

	for (int i = 0; i < n; i++) {
		otherArray[i] = i; 
	}
	double max = maxError(myArray, otherArray, n);
	double expectedMax = std::log10(std::fabs(n-1 - 1));
	if (fabs(expectedMax - max) > eps) {
		success = false;
	}

	otherArray[10] = 50;
	max = maxError(myArray, otherArray, n);
	expectedMax = std::log10(std::fabs(50 - 1));
	if (fabs(expectedMax - max) > eps) {
		success = false;
	}
	if (success) {
		printf("\t- PASSED : max value matched expected max value\n\n");
	} else {
		printf("\t- FAILED : max value did not match expected max value\n\n");
	}

	printf("\nvoid generalTriDiaSolver(double *a, double *b, double *c, double *v, double *f, int n)\n");

	n = 20;
	double *a, *b, *c, *x, *v, *f; 
	srand(time(NULL));
	v = createVector(0, n);
	f = createVector(0, n);
	x = createVector(0, n);
	a = createVector(0, n-1);
	b = createVector(0, n);
	c = createVector(0, n-1);
	float multiplier = 1;

	for (int i = 0; i < (n - 1); i++) {
		x[i] = rand()/double(RAND_MAX)*multiplier;
		a[i] = rand()/double(RAND_MAX)*multiplier;
		b[i] = rand()/double(RAND_MAX)*multiplier;
		c[i] = rand()/double(RAND_MAX)*multiplier;

	}
	x[n-1] = rand()/double(RAND_MAX)*multiplier;
	b[n-1] = rand()/double(RAND_MAX)*multiplier;

	f[0] = x[0]*b[0] + x[1]*c[0];
	for (int i = 1; i < (n - 1); i++) {
		f[i] = a[i-1]*x[i-1] + b[i]*x[i] + c[i]*x[i+1];
	}
	f[n-1] = a[n-2]*x[n-2] + b[n-1]*x[n-1];

	generalTriDiaSolver(a, b, c, v, f, n);

	success = true;
	eps = 1e-12;

	for (int i = 0; i < n; i++) {
		if (fabs(v[i] - x[i]) > eps) {
			success = false;
			break;
		}
	}
	if (success) {
		printf("\t- PASSED : calculated array v seems correct\n\n");
	} else {
		printf("\t- FAILED : calculated array v does not match expected\n\n");
	}

	for (int i = 0; i < n; i++) {printf("%.30f\n", fabs(v[i] - x[i]));}


	printf("\nvoid constTriDiaSolver(double a, double c, double *b, double *v, double *f, int n)\n");


	return 0;
}













