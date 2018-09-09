#include <cmath>
#include "armadillo"
#include <stdlib.h>
#include <string>

using namespace std;
using namespace arma;

double *createVector(double value, int n) {
	//  Allocating space for a vector and fill elements with a constant value

	double *vector;
	vector = new double[n];
	for (int i = 0; i < n; i++) {
		vector[i] = value;
	}
	return vector;
}

double *solutionVector(int n) {
	// Function for creating the right hand side vector of the linear equation 

	double *f;
	f = new double[n];

	double h = 1.0/(n + 1.0);
	double hh = h*h;
	double x = h;
	for (int i = 0; i < n; i++) {
		x += h;
		f[i] = 100.0*exp(-10.0*x)*hh;
	}
	return f;
}

void generalTriDiaSolver(double *a, double *b, double *c, double *v, double *f, int n) {
	// General algorithm for performing TDMA 

	double *d, *g;
	double temp;
	d = createVector(b[0], n);
	g = createVector(f[0], n);
	// Forward substitiution
	for (int i = 1; i < n; i++) {
		temp = a[i]/d[i-1];
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

void constTriDiaSolver(double a, double c, double *b, double *v, double *f, int n) {
	//Special algorithm with constant values

	double *d, *g;
	d = createVector(b[0], n);
	g = createVector(f[0], n);
	double ac = a*c;
	// Forward substitution
	for (int i = 1; i < n; i++) {
		d[i] = b[i] - ac/d[i-1];
		g[i] = f[i] - a*g[i-1]/d[i-1];
	}
	//Backward substitution
	v[n-1] = g[n-1]/d[n-1];
	for (int j = n - 2; j > -1; j--) {
		v[j] = (g[j] - c*v[j+1])/d[j];
	}
	delete[] d;
	delete[] g;
}

double *exactSolution(int n) {
	// Create vector with analytical solution

	double *u;
	u = new double[n];

	double h = 1.0/(n + 1.0);
	double x = h;
	for (int i = 0; i < n; i++) {
		x += h;
		u[i] = 1 - (1 - exp(-10.0))*x - exp(-10.0*x);
	}
	return u;
}

vec armadillo_LU_solve(mat A, vec b) {
	mat L, U;
	lu(L, U, A);

	// only for optimization 
	mat L_prime = trimatl(L);
	mat U_prime = trimatu(U);

	//To solve LUx = b, you first solve Ly = b and then Ux = y.
	vec y, x;
	y = solve(L_prime, b);
	x = solve(U_prime, y);

	return x;
}	


