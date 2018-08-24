#include "pch.h"
#include <iostream>
#include <cmath>

double * createVector(int n, double value);
void generalTriDiaSolver(double * a, double * b, double * c, double * v, double * f, int n);
void constTriDiaSolver(int a, int c, double * b, double * v, double * f, int n);
double * exactSolution(int n);

int main(int argc, char* argv[])
{
	// Reading in dimensons and values for vectors
	int n = atoi(argv[1]);
	double a_value = atoi(argv[2]);
	double b_value = atoi(argv[3]);
	double c_value = atoi(argv[4]);

	// Create vectors
	double *a, *b, *c, *f, *v, *exact;
	a = createVector(n+2, a_value);
	b = createVector(n+2, b_value);
	c = createVector(n+2, c_value);
	f = createVector(n+2, 0.0);
	v = createVector(n+2, 0.0);
	exact = exactSolution(n);

	// Fill the solution vector
	double h = 1.0/(n + 1);
	for (int i = 1; i < n+1; i++) {
		f[i] = 100*std::exp(-10.0*i*h)*pow(h, 2);
	}

	generalTriDiaSolver(a, b, c, v, f, n);
	//constTriDiaSolver(a_value, c_value, b, u, f, n);

	double error;
	for (int i = 0; i < n+2; i++) {
		error = v[i] - exact[i];
		printf("v: %.5f,   u: %.10f,   %f\n", exact[i], v[i], error);
	}
	printf("last exact: %f", exact[n-1]);


	delete[] a;
	delete[] b;
	delete[] c;
	delete[] f;
	delete[] v;
	return 0;
}

//  Allocating space for a vector and fill elements with a constant value
double * createVector(int n, double value) {
	double * vec;
	vec = new double[n];
	for (int i = 0; i < n; i++) {
		vec[i] = value;
	}
	return vec;
}

// General algorithm for performing TDMA 
void generalTriDiaSolver(double * a, double * b, double * c, double * v, double * f, int n) {
	// Forward substitiution
	for (int i = 1; i < (n+2); i++) {
		b[i] = b[i] - a[i]*c[i-1]/b[i-1];
		f[i] = f[i] - a[i]*f[i-1]/b[i-1];
	}

	v[n] = f[n] / b[n];
	// Backward substitution 
	for (int j = n; j > 0; j--) {
		v[j] = (f[j] - c[j]*v[j+1])/b[j];
	}
}

//Special algorithm with constant values
void constTriDiaSolver(int a, int c, double * b, double * v, double * f, int n) {
	double m = a*a;
	// Forward substitution
	for (int i = 1; i < n; i++) {
		b[i] = b[i] - m/b[i-1];
		f[i] = f[i] - a*f[i-1]/b[i-1];
	}
	v[n] = v[n]/b[n];
	//Backward substitution
	for (int j = n-1; j > 0; j--) {
		v[j] = (f[j] - c*v[j+1])/b[j];
	}
}
// Create vector with analytical solution
double * exactSolution(int n) {
	double * u;
	double x;
	u = new double[n];
	double h = 1.0/(n + 1);
	for (int i = 0; i < n+1; i++) {
		x = i*h;
		u[i] = 1 - (1 - std::exp(-10))*x - std::exp(-10 * x);
	}
	return u;
}


