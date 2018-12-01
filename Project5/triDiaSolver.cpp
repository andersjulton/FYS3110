#include "utils.h"

void generalTriDiaSolver(double *a, double *b, double *c, double *v, double *f, int n) {
	// General algorithm for performing TDMA 

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

void constTriDiaSolver(double a, double b, double c, double *v, double *f, int n) {
	//Special algorithm with constant values
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
