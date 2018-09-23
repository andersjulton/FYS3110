#include <fstream>
#include "func.h"
#include "utils.h"
#include "armadillo"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <string>

using namespace std;
using namespace arma;

// Solving Ax = b
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

//Testing the functions for Thomas algorithm
double maxErrorDiaSolver(int n) {

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
	error = maxEpsilon(u, v, n);

	delete[] f;
	delete[] v;
	delete[] u;
	return error;
}

//Compare TriDiaSolver general and with constants
void compareTime(int a_value, int b_value, int c_value, int n) {

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

//Compare TriDiaSolver with Armadillo
void compareTimeArmadillo(int n) {
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