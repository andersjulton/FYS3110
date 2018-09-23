#include <cmath>
#include <stdlib.h>
#include "utils.h"
#include <string>
#include <algorithm>
#include <fstream>

using namespace std;


// Function for creating the right hand side vector of the linear equation 
double *solutionVector(int n) {
	double *f;
	f = new double[n];

	double h = 1.0/(n + 1.0);
	double hh = h*h;
	double x = h;
	for (int i = 0; i < n; i++) {
		f[i] = 100.0*exp(-10.0*x)*hh;
		x += h;
	}
	return f;
}

// General algorithm for performing TDMA 
void generalTriDiaSolver(double *a, double *b, double *c, double *v, double *f, int n) {
	double *d, *g;
	double temp;
	d = createVector(b[0], n);
	g = createVector(f[0], n);
	// Forward substitiution
	for (int i = 1; i < n; i++) {
		temp = a[i-1]/d[i-1];
		d[i] = b[i] - temp*c[i-1];
		g[i] = f[i] - temp*g[i-1];
	}
	// Backward substitution 
	v[n-1] = g[n-1]/d[n-1];
	for (int j = n-2; j > -1; j--) {
		v[j] = (g[j] - c[j]*v[j+1])/d[j];
	}
	delete[] d;
	delete[] g;
}

//Special algorithm with constant values
void constTriDiaSolver(double a, double b, double c, double *v, double *f, int n) {
	double *d, *g;
	d = createVector(b, n);
	g = createVector(f[0], n);
	double ac = a*c;
	// Forward substitution
	for (int i = 1; i < n; i++) {
		d[i] = b - ac/d[i-1];
		g[i] = f[i] - g[i-1]*a/d[i-1];
	}
	//Backward substitution
	v[n-1] = g[n-1]/d[n-1];
	for (int j = n - 2; j > -1; j--) {
		v[j] = (g[j] - c*v[j+1])/d[j];
	}
	delete[] d;
	delete[] g;
}

// Create vector with analytical solution
double *exactSolution(int n) {
	double *u;
	u = new double[n];

	double h = 1.0/(n + 1.0);
	double x = h;
	for (int i = 0; i < n; i++) {
		u[i] = 1.0 - (1.0 - exp(-10.0))*x - exp(-10.0*x);
		x += h;
	}
	return u;
}
