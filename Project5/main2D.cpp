#define _USE_MATH_DEFINES
#include "PDE.h"
#include "utils.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void writeToFile(double t, double h, string filename);

int main() {
	double t1 = 0.05;
	double t2 = 0.005;
	printf("t = %.2f, t = %.2f\n", t1, t2);
	writeToFile(t1, 0.01, "t1_h_0.01");
	writeToFile(t2, 0.01, "t2_h_0.01");
	writeToFile(t1, 0.1, "t1_h_0.1");
	writeToFile(t2, 0.1, "t2_h_0.1");
	return 0;
}

void writeToFile(double t, double h, string filename){
	double dt = 0.25*(h*h*h*h);					// h = dx = dy
	printf("h = %f, t = %f \n", h, t);
	int timeSteps = (int) (t/dt + 1);
	int n = (int) (1.0/h + 1);				// L = 1
	int m = n;
	printf("n = %d\n", n);
	double alpha = dt/(h*h);
	printf("time steps = %d\n", timeSteps);

	double **u = createMatrix(n, n);
	double **exact = createMatrix(n, n);

	double x, y;
	// edges filled with zero
	for (int i = 1; i < (m-1); i++) {
		y = h*i;
		for (int j = 1; j < (n-1); j++) {
			x = h*j;
			// using-ish eq. from page 314, just needed something to test, no idea if it's correct
			// using L = 1 & n = 1
			u[i][j] = sin(M_PI*x)*sin(M_PI*y);						// ------FIX THIS-------
			exact[i][j] = u[i][j]*exp(-2*M_PI*M_PI*t);				// ------FIX THIS-------
		};
	}
	forwardEuler(u, alpha, timeSteps, m, n);
	double **error = absError(exact, u, n, n);

	doubleMatrixToFile(u, n, n, filename);
	doubleMatrixToFile(error , n, n, filename + "_error");

	deleteMatrix(u, n);
	deleteMatrix(exact, n);
	deleteMatrix(error, n);;
}
