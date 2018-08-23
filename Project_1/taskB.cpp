#include <math.h>
#include <stdio.h>

void triaSolve(double *a, double *b, double *c, double *v, double *f, int n) {
	// with (way too) much help from "Teach yourself C++". 
	// not even sure if this is correct
	for (int i = 1; i < (n+2); i++) {
		//printf("b: %f ", b[i]);
		b[i] = b[i] - c[i-1]*a[i-1]/b[i-1];
		//printf("b: %f ", b[i]);
		f[i] = f[i] - f[i-1]*a[i-1]/b[i-1]; // ????
		//printf("f: %f ", f[i]);
	}
	//v[n-1] = f[n-1]/b[n-1];
	printf("\n ");

	for (int i = n; i > 0; i--) {
		f[i] = f[i] - f[i+1]*c[i]/b[i];
		//printf("f: %f ", f[i]);
		v[i] = f[i]/b[i];
		//printf("v: %f ", v[i]);
	}
	printf("\n ");
}

double * makeVector( int n, double value) {
	double *vector;
	vector = new double[n];
	for (int i = 0; i < n; i++) {
		vector[i] = value;
	}
	return vector;
}

void exact(double *u, double *f, int n) {
	double x, h;

	h = 1.0/(n + 1);		// I hate this, I want to change it to h/n, and just index differently
	printf("%f", h);

	for (int i = 1; i < n+1; i++) {
		x = i*h;
		double z = 1.0 - (1 - exp(-10))*x - exp(-10*x);
		//printf("x is %f, z is %.40f", x, z);
		u[i] = z;
		f[i] = 100*exp(-10.0*x)*h*h;
	}
}

int main() {
	// just checking if this even works...

	double *a, *b, *c, *f, *v, *u;
	int n = 100;

	a = makeVector(n+2, -1.0);
	b = makeVector(n+2, 2.0);
	c = makeVector(n+2, -1.0);
	f = makeVector(n+2, 0);
	v = makeVector(n+2, 0);
	u = makeVector(n+2, 0);

	exact(u, f, n);
	triaSolve(a, b, c, v, f, n);

	double error;
	for (int i = 0; i < n+2; i++) {
		error = u[i] - v[i];
		printf("v: %.5f,   u: %.40f,   %f\n", v[i], u[i], error);
	}

	delete[] u;
	delete[] a;
	delete[] b;
	delete[] c;
	delete[] f;
	delete[] v;
	return 0;
}

/*
Eirills-MacBook-Pro:Project_1 eirillsh$ c++ -c -Wall taskB.cpp && c++ -o taskB.exe  taskB.o && ./taskB.exe
Not sure if this is correct in any way. The results are pretty off. They get better as n gets higher.

*/


