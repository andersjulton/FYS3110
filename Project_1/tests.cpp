#include "thomas.h"
#include <cmath>
#include <chrono>
#include <algorithm>
#include "armadillo"

using namespace arma;
using namespace std;

double maxError(double *u, double *v, int n) {
	// Find the largest relative error in array v in respect to array u

	double maxError = std::log10(std::fabs((v[0] - u[0])/u[0]));
	double epsilon;
	for (int i = 1; i < n; i++) {
		epsilon = std::log10(std::fabs((v[i] - u[i])/u[i]));
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
	double *b, *f, *v, *u;
	b = createVector(b_value, n);
	v = createVector(0.0, n);

	u = exactSolution(n);
	f = solutionVector(n);

	constTriDiaSolver(a_value, c_value, b, v, f, n);

	error = maxError(u, v, n);

	delete[] b;
	delete[] f;
	delete[] v;
	delete[] u;
	return error;
}

void compareTime(int a_value, int b_value, int c_value, int n) {
	//Compare TriDiaSolver general and with constants

	std::chrono::duration<double> elapsed;
	double *timeGeneralTriDiaSlover, *timeConstTriDiaSlover;

	timeGeneralTriDiaSlover = new double[5];
	timeConstTriDiaSlover = new double[5];

	double *a, *b, *c, *f, *v;
	a = createVector(a_value, n);
	b = createVector(b_value, n);
	c = createVector(c_value, n);
	v = createVector(0.0, n);
	f = solutionVector(n);

	printf("For n = %d\n", n);
	for (int i = 0; i < 5; i++) {

		auto begin = std::chrono::high_resolution_clock::now();
		generalTriDiaSolver(a, b, c, v, f, n);
		auto end = std::chrono::high_resolution_clock::now();
		elapsed = (end - begin);
		timeGeneralTriDiaSlover[i] = (double)elapsed.count();

		begin = std::chrono::high_resolution_clock::now();
		constTriDiaSolver(a_value, c_value, b, v, f, n);
		end = std::chrono::high_resolution_clock::now();
		elapsed = (end - begin);
		timeConstTriDiaSlover[i] = (double)elapsed.count();
	}

	std::sort(timeGeneralTriDiaSlover, timeGeneralTriDiaSlover + 5);
	printf("Time it took for generalTriDiaSolver: %g s\n", timeGeneralTriDiaSlover[2]);

	std::sort(timeConstTriDiaSlover, timeConstTriDiaSlover + 5);
	printf("Time it took for constTriDiaSolver: %g s\n", timeConstTriDiaSlover[2]);

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

	std::chrono::duration<double> elapsed;
	double *timeTriDiaSlover, *timeArmadillo;

	timeTriDiaSlover = new double[5];
	timeArmadillo = new double[5];

	double a_value = -1.0;
	double b_value = 2.0;
	double c_value = -1.0;

	//Create vectors
	double *b, *f, *v;
	b = createVector(b_value, n);
	v = createVector(0.0, n);
	f = solutionVector(n);

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

	std::vector<double> f_vector;
	f_vector.assign(f, f + n);
	vec y(f_vector);

	printf("\nFor n = %d\n", n);
	for (int i = 0; i < 5; i++) {
		// own code
		auto begin = std::chrono::high_resolution_clock::now();
		constTriDiaSolver(a_value, c_value, b, v, f, n);
		auto end = std::chrono::high_resolution_clock::now();
		elapsed = (end - begin);
		timeTriDiaSlover[i] = (double)elapsed.count();
		printf("Time elapsed: %.7f seconds\n", timeTriDiaSlover[i]);
		// armadillo
		begin = std::chrono::high_resolution_clock::now();
		x = armadillo_LU_solve(A, y);
		end = std::chrono::high_resolution_clock::now();
		elapsed = (end - begin);
		timeArmadillo[i] = (double)elapsed.count();
		printf("Time elapsed: %.7f seconds\n", timeArmadillo[i]);

		double maxError = fabs(v[0] - x(0))/x(0);
		double error;
		for (int i = 1; i < n; i++) {
			error = fabs(v[i] - x(i))/x(i);
			if (error > maxError) {
				maxError = error;
			}
		}
		printf("Largest relative error was %.20f\n", (maxError*100));
	}

	std::sort(timeTriDiaSlover, timeTriDiaSlover + 5);
	printf("Time it took for TriDiaSolver: %g s\n", timeTriDiaSlover[2]);

	std::sort(timeArmadillo, timeArmadillo + 5);
	printf("Time it took for armadillo: %g s\n", timeArmadillo[2]);

	delete[] b;
	delete[] f;
	delete[] v;
	delete[] timeTriDiaSlover;
	delete[] timeArmadillo; 
}
