#include <cmath>
#include "armadillo"
#include <stdlib.h>
#include <string>
#include <chrono>
#include <algorithm>
#include <fstream>

using namespace std;
using namespace arma;

void printError(double *u, double *v, int n) {
	// write out difference between u-array and v-array

	double error;
	printf("\n");
	for (int i = 0; i < n; i++) {
		error = u[i] - v[i];
		printf("v = %.20f, u = %.20f, error = %.30f\n", v[i], u[i], error);
	}
	printf("\n");
}

void writeToFile(double *v, int n) {
	// write v-array to a txt-file

	ofstream myfile("n_"+to_string(n)+".txt");
	if (myfile.is_open()) {
		myfile << (n + 2) << "\n";
		myfile << 0 << "\n";
		for (int i = 0; i < n; i++) {
			myfile << v[i] << "\n";
		}
		myfile << 0 << "\n";
	}
}

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
		f[i] = 100.0*exp(-10.0*x)*hh;
		x += h;
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

double *exactSolution(int n) {
	// Create vector with analytical solution

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

vec armadillo_LU_solve(mat A, vec b) {
	// Solving Ax = b

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

double maxError(double *u, double *v, int n) {
	// Find the largest relative error in array v in respect to array u

	double maxError = log10(fabs((v[0] - u[0])/u[0]));
	double epsilon;
	for (int i = 1; i < n; i++) {
		epsilon = log10(fabs((v[i] - u[i])/u[i]));
		if (maxError < epsilon) {
			maxError = epsilon;
		}
	}
	return maxError;
}

double maxErrorDiaSolver(int n) {
	//Testing the functions for Thomas algorithm

	double a_value = -1.0;
	double b_value = 2.0;
	double c_value = -1.0;
	double error;

	// Create vectors
	double *f, *v, *u;
	v = createVector(0.0, n);
	u = exactSolution(n);
	f = solutionVector(n);

	constTriDiaSolver(a_value, b_value, c_value, v, f, n);
	error = maxError(u, v, n);

	delete[] f;
	delete[] v;
	delete[] u;
	return error;
}

void compareTime(int a_value, int b_value, int c_value, int n) {
	//Compare TriDiaSolver general and with constants

	chrono::duration<double> elapsed;
	double *timeGeneralTriDiaSlover, *timeConstTriDiaSlover;

	timeGeneralTriDiaSlover = new double[5];
	timeConstTriDiaSlover = new double[5];

	double *a, *b, *c, *f, *v;
	a = createVector(a_value, n-1);
	b = createVector(b_value, n);
	c = createVector(c_value, n-1);
	v = createVector(0.0, n);
	f = solutionVector(n);

	printf("For n = %d\n", n);
	// timing 5 times
	for (int i = 0; i < 5; i++) {

		auto begin = chrono::high_resolution_clock::now();
		generalTriDiaSolver(a, b, c, v, f, n);
		auto end = chrono::high_resolution_clock::now();
		elapsed = (end - begin);
		timeGeneralTriDiaSlover[i] = (double)elapsed.count();

		begin = chrono::high_resolution_clock::now();
		constTriDiaSolver(a_value, b_value, c_value, v, f, n);
		end = chrono::high_resolution_clock::now();
		elapsed = (end - begin);
		timeConstTriDiaSlover[i] = (double)elapsed.count();
	}
	// Finding the median
	sort(timeGeneralTriDiaSlover, timeGeneralTriDiaSlover + 5);
	printf("Time it took for generalTriDiaSolver: %g s\n", timeGeneralTriDiaSlover[2]);

	sort(timeConstTriDiaSlover, timeConstTriDiaSlover + 5);
	printf("Time it took for constTriDiaSolver: %g s\n", timeConstTriDiaSlover[2]);

	printf("Ratio: %g \n", (timeGeneralTriDiaSlover[2]/timeConstTriDiaSlover[2]));

	delete[] a;
	delete[] b;
	delete[] c;
	delete[] f;
	delete[] v;
	delete[] timeGeneralTriDiaSlover;
	delete[] timeConstTriDiaSlover;
}

void compareTimeArmadillo(int n) {
	//Compare TriDiaSolver with Armadillo

	chrono::duration<double> elapsed;
	double *timeTriDiaSlover, *timeArmadillo;

	timeTriDiaSlover = new double[5];
	timeArmadillo = new double[5];

	double a_value = -1.0;
	double b_value = 2.0;
	double c_value = -1.0;

	//Create vectors
	double *f, *v;
	v = createVector(0.0, n);
	f = solutionVector(n);

	// create A in Ax = y
	mat A(n, n, fill::zeros);
	A(0, 0) = b_value;
	A(0, 1) = c_value;
	for (int i = 1; i < (n - 1); i++) {
		A(i, i-1) = a_value;
		A(i, i) = b_value;
		A(i, i+1) = c_value;
	}
	A(n-1, n-2) = a_value;
	A(n-1, n-1) = b_value;

	vec x;
	vector<double> f_vector;
	f_vector.assign(f, f + n);
	vec y(f_vector);

	printf("\nFor n = %d\n", n);
	// timing 5 times
	for (int i = 0; i < 5; i++) {
		// time own code
		auto begin = chrono::high_resolution_clock::now();
		constTriDiaSolver(a_value, b_value, c_value, v, f, n);
		auto end = chrono::high_resolution_clock::now();
		elapsed = (end - begin);
		timeTriDiaSlover[i] = (double)elapsed.count();
		// time armadillo
		begin = chrono::high_resolution_clock::now();
		x = armadillo_LU_solve(A, y);
		end = chrono::high_resolution_clock::now();
		elapsed = (end - begin);
		timeArmadillo[i] = (double)elapsed.count();
	}
	// Finding the median
	sort(timeTriDiaSlover, timeTriDiaSlover + 5);
	printf("Time it took for TriDiaSolver: %g s\n", timeTriDiaSlover[2]);

	sort(timeArmadillo, timeArmadillo + 5);
	printf("Time it took for armadillo: %g s\n", timeArmadillo[2]);

	printf("Ratio: %g \n", timeArmadillo[2]/timeTriDiaSlover[2]);

	delete[] f;
	delete[] v;
	delete[] timeTriDiaSlover;
	delete[] timeArmadillo; 
}	
