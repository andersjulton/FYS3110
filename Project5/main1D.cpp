#include "PDE.h"
#include "utils.h"
#include "exact.h"
#include <iostream>
#include <fstream>

using namespace std;

void writeToFile(double t, double dx, string filename);

void stability(int calc, int start, double incr, string name);

void dx(int calc, int start, double incr);

int main() {
	double t1 = 0.05;
	double t2 = 0.5;
	printf("t = %.2f, t = %.2f\n", t1, t2);
	writeToFile(t1, 0.01, "t1_dx_0.01");
	writeToFile(t2, 0.01, "t2_dx_0.01");
	writeToFile(t1, 0.1, "t1_dx_0.1");
	writeToFile(t2, 0.1, "t2_dx_0.1");

	stability(5, 1, 0.1, "1");
	stability(5, 10, 0.05, "2");

	dx(1000, 20, 0.001);
	return 0;
}

void writeToFile(double t, double dx, string filename) {
	double dt = 0.5*(dx*dx);
	int timeSteps = (int) (t/dt + 1);
	int n = (int) (1.0/dx + 1);				// L = 1
	double alpha = dt/(dx*dx);

	double **u = createMatrix(3, n);
	double *exact = createVector(1, n); 
	double **error = new double*[3];

	analytic1D(exact, t, dx, n);

	u[0][n-1] = 1;
	u[1][n-1] = 1;
	u[2][n-1] = 1;

	forwardEuler(u[0], alpha, timeSteps, n);
	backwardEuler(u[1], alpha, timeSteps, n);
	crank_nicolson(u[2], alpha, timeSteps, n);

	error[0] = absError(exact, u[0], n);
	error[1] = absError(exact, u[1], n);
	error[2] = absError(exact, u[2], n);


	ofstream outfile, errorfile;
	outfile.open(filename + ".txt");
	errorfile.open(filename + "_error.txt");

	outfile << "Forward Euler\tBackward Euler\tCrank-Nicolson\texact";
	errorfile << "Forward Euler\tBackward Euler\tCrank-Nicolson";

	for (int i = 0; i < n; i++) {
		outfile << "\n";
		errorfile << "\n";
		for (int j = 0; j < 3; j++) {
			outfile << to_string(u[j][i]) + "\t\t ";
			errorfile << to_string(error[j][i]) + "\t\t ";
		}
		outfile << to_string(exact[i]);
	}
	outfile.close();
	errorfile.close();

	delete[] exact;
	deleteMatrix(u, 3);
	deleteMatrix(error, 3);
}

void stability(int calc, int start, double incr, string name) {
	double t = 0.05;
	double dx = 0.1;
	double dt, fac;
	int timeSteps;
	int n = (int) (1.0/dx + 1);				// L = 1
	printf("%d\n", n);
	double alpha;

	double **u = createMatrix(calc, n);
	double *exact = createVector(1, n);
	double **error = new double*[calc];

	analytic1D(exact, t, dx, n);

	string filename = "stability" + name;
	ofstream outfile, errorfile;
	outfile.open(filename + ".txt");
	errorfile.open(filename + "_error.txt");

	for (int i = 0; i < calc; i++) {
		u[i][n-1] = 1;
		fac = incr*(i + start);
		dt = fac*(dx*dx);
		alpha = dt/(dx*dx);
		timeSteps = (int) (t/dt + 1);

		outfile << to_string(fac) + "\t";
		errorfile << to_string(fac) + "\t";

		forwardEuler(u[i], alpha, timeSteps, n);
		error[i] = absError(exact, u[i], n);
	}
	outfile << "exact";

	for (int i = 0; i < n; i++) {
		outfile << "\n";
		errorfile << "\n";
		for (int j = 0; j < calc; j++) {
			outfile << to_string(u[j][i]) + "\t ";
			errorfile << to_string(error[j][i]) + "\t ";
		}
		outfile << to_string(exact[i]);
	}
	outfile.close();
	errorfile.close();

	delete[] exact;
	deleteMatrix(u, calc);
	deleteMatrix(error, calc);
}

void dx(int calc, int start, double incr) {
	double t = 0.05;
	double dx = 0.1;
	double dt, fac;
	int timeSteps;
	int n = (int) (1.0/dx + 1);				// L = 1
	printf("%d\n", n);
	double alpha;

	double *u;
	double *exact = createVector(1, n);
	double error;

	analytic1D(exact, t, dx, n);

	string filename = "dt_error.txt";
	ofstream outfile;
	outfile.open(filename);

	for (int i = 0; i < calc; i++) {
		u = createVector(0, n);
		u[n-1] = 1;

		fac = incr*(i + start);
		dt = fac*(dx*dx);
		alpha = dt/(dx*dx);
		timeSteps = (int) (t/dt + 1);
		outfile << to_string(fac) + "\t";

		forwardEuler(u, alpha, timeSteps, n);
		error = maxAbsError(exact, u, n);
		outfile << to_string(error) + "\n";

		delete[] u;
	}
	outfile.close();

	delete[] exact;
}