#include "tests.h"
#include <armadillo>
#include "thomas.h"

// Find the largest relative error in array v in respect to array u
double maxError(double *u, double *v, int n) {
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

//Testing the functions for Thomas algorithm
double maxErrorDiaSolver(int n) {
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

//Compare TriDiaSolver general and with constants
void compareTime(int a_value, int b_value, int c_value, int n) {
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
		timeGeneralTriDiaSlover[i] = (double) elapsed.count();
		
		begin = std::chrono::high_resolution_clock::now();
		constTriDiaSolver(a_value, c_value, b, v, f, n);
		end = std::chrono::high_resolution_clock::now();
		elapsed = (end - begin);
		timeConstTriDiaSlover[i] = (double) elapsed.count();
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

//Compare TriDiaSolver with Armadillo
void compareTimeArmadillo(int n) {
	std::chrono::duration<double> elapsed;
	double *timeTriDiaSlover, *timeArmadillo;

	timeTriDiaSlover = new double[5];
	timeArmadillo = new double[5];

	double a_value = -1.0;
	double b_value = 2.0;
	double c_value = -1.0;

	// Create vectors
	double *b, *f, *v;
	b = createVector(b_value, n);
	v = createVector(0.0, n); 
	f = solutionVector(n);

	printf("\nFor n = %d\n", n);
	for (int i = 0; i < 5; i++) {

		auto begin = std::chrono::high_resolution_clock::now();
		constTriDiaSolver(a_value, c_value, b, v, f, n);
		auto end = std::chrono::high_resolution_clock::now();
		elapsed = (end - begin);
		//std::cout << "Time elapsed: " << elapsed.count() << " seconds" << std::endl;
		timeTriDiaSlover[i] = (double) elapsed.count();
		
		begin = std::chrono::high_resolution_clock::now();
		// !!! What?
		end = std::chrono::high_resolution_clock::now();
		elapsed = (end - begin);
		//std::cout << "Time elapsed: " << elapsed.count() << " seconds" << std::endl;
		timeArmadillo[i] = (double) elapsed.count();
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
