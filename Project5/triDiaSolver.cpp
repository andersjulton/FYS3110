#include "utils.h"

//Special algorithm with constant values with boundary conditions
void triDiaSolver(double a, double b, double c, double *v, double *f, int n) {
	double *d, *g;
	d = createVector(b, n);
	g = createVector(f[0], n);
	double ac = a*c;
	// Forward substitution
	for (int i = 1; i < (n-1); i++) {
		d[i] = b - ac/d[i-1];
		g[i] = f[i] - g[i-1]*a/d[i-1];
	}
	//Backward substitution
	v[n-1] = f[n-1];
	for (int j = (n-2); j > 0; j--) {
		v[j] = (g[j] - c*v[j+1])/d[j];
	}
	v[0] = f[0];
	delete[] d;
	delete[] g;
}