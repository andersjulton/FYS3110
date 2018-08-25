//#include "pch.h"    do we need this?
#include <iostream>
#include <cmath>

double * createVector(double value, int n);
double * solutionVector(int n);
void generalTriDiaSolver(double * a, double * b, double * c, double * v, double * f, int n);
void constTriDiaSolver(int a, int c, double * b, double * v, double * f, int n);
double * exactSolution(int n);
double maxErrorDiaSolver(double *u, double *v, int n);
void testDiaSolver(int n);
void printError(double *u, double *v, int n);

int main(int argc, char* argv[]) {
	// Reading in dimensons and values for vectors
	int n = atoi(argv[1]);
	double a_value = atoi(argv[2]);
	double b_value = atoi(argv[3]);
	double c_value = atoi(argv[4]);

	// Create vectors
	double *b, *f, *v, *u;
	double *a, *c;

	a = createVector(a_value, n);
	b = createVector(b_value, n);
	c = createVector(c_value, n);
	v = createVector(0.0, n);

	u = exactSolution(n);   // this is only an exact solution for ONE case
	f = solutionVector(n);

	generalTriDiaSolver(a, b, c, v, f, n);
	//constTriDiaSolver(a_value, c_value, b, u, f, n);

	printError(u, v , n);

	delete[] a;
	delete[] b;
	delete[] c;
	delete[] f;
	delete[] v;
	delete[] u;


	//testing:
	bool testing = true;
	if (testing) {
		testDiaSolver(10);
		testDiaSolver(100);
		testDiaSolver(1000);
		testDiaSolver(10000000);
	}
	return 0;
}

//  Allocating space for a vector and fill elements with a constant value
double * createVector(double value, int n) {
	double *vector;
	vector = new double[n];
	for (int i = 0; i < n; i++) {
		vector[i] = value;
	}
	return vector;
}

double * solutionVector(int n) {
	// Fill the solution vector
	double *f;
	f = new double[n];
	double h = 1.0/(n + 1.0);
	for (int i = 1; i < n+1; i++) {
		f[i-1] = 100.0*std::exp(-10.0*i*h)*h*h; 
	}
	return f;
}

// General algorithm for performing TDMA 
void generalTriDiaSolver(double *a, double *b, double *c, double *v, double *f, int n) {
	// Forward substitiution
	for (int i = 1; i < n; i++) {
		b[i] = b[i] - a[i]*c[i-1]/b[i-1]; 
		f[i] = f[i] - a[i]*f[i-1]/b[i-1];  
	}
	// Backward substitution 
	v[n-1] = f[n-1]/b[n-1];
	for (int j = n-2; j > -1; j--) {
		v[j] = (f[j] - c[j]*v[j+1])/b[j];
	}
}

//Special algorithm with constant values
void constTriDiaSolver(int a, int c, double *b, double *v, double *f, int n) {
	double ac = a*c;
	// v[0] = v[n+1] is already 0.
	// Forward substitution
	for (int i = 1; i < n; i++) {
		b[i] = b[i] - ac/b[i-1];
		f[i] = f[i] - a*f[i-1]/b[i-1];
	}
	//Backward substitution
	v[n-1] = f[n-1]/b[n-1];
	for (int j = n - 2; j > -1; j--) {
		v[j] = (f[j] - c*v[j+1])/b[j];
	}
}

// Create vector with analytical solution
double * exactSolution(int n) {
	double x;
	double h = 1.0/(n + 1.0);
	double *u;
	u = new double[n];

	for (int i = 1; i < n+1; i++) { 
		x = i*h;
		u[i-1] = 1 - (1 - std::exp(-10.0))*x - std::exp(-10.0*x);
	}
	return u;
}

double maxErrorDiaSolver(double *u, double *v, int n) {
	double maxError = std::log10( std::fabs( (v[1] - u[1])/u[1] ));
	//double maxError = abs( (v[1] - u[1])/u[1] );
	double epsilon;

	for (int i = 2; i < n+1; i++) {
		epsilon = std::log10( std::fabs( (v[i] - u[i])/u[i] ));
		//epsilon = abs( (v[i] - u[i])/u[i] );
		if (maxError < epsilon) {
			maxError = epsilon;
		}
	}
	return maxError;
}

void testDiaSolver(int n) {
	double a_value = -1.0;
	double b_value = 2.0;
	double c_value = -1.0;
	double maxError;

	// Create vectors
	double *b, *f, *v, *u;
	double *a, *c;
	a = createVector(a_value, n);
	b = createVector(b_value, n);
	c = createVector(c_value, n);
	v = createVector(0.0, n);

	u = exactSolution(n);   // this is only an exact solution for ONE case
	f = solutionVector(n);

	//generalTriDiaSolver(a, b, c, v, f, n);
	constTriDiaSolver(a_value, c_value, b, v, f, n);

	maxError = maxErrorDiaSolver(u, v, n);
	printf("For n = %10d: the max value of relative error is %.20f\n", n, maxError);

	delete[] a;
	delete[] b;
	delete[] c;
	delete[] f;
	delete[] v;
	delete[] u;
}

void printError(double *u, double *v, int n) {
	// for debugging
	double error;
	printf("\n");
	for (int i = 0; i < n; i++) {
		error = u[i] - v[i];
		printf("v = %.20f, u = %.20f, error = %.30f\n", v[i], u[i], error);
	}
	printf("\n");
}

// eirillsh$ c++ -c -Wall Project_1.cpp && c++ -o Project_1.exe  Project_1.o && ./Project_1.exe 10 -1 2 -1



