#include "PDE.h"
#include "utils.h"
#include <iostream>
#include <fstream>

using namespace std;

void writeToFile(double t, double dx, string filename);

int main() {
	double t1 = 0.5;
	double t2 = 0.05;
	printf("t = %.2f, t = %.2f\n", t1, t2);
	writeToFile(t1, 0.01, "t1_dx_0.01");
	writeToFile(t2, 0.01, "t2_dx_0.01");
	writeToFile(t1, 0.1, "t1_dx_0.1");
	writeToFile(t2, 0.1, "t2_dx_0.1");
	return 0;
}

void writeToFile(double t, double dx, string filename) {
	double dt = 0.5*(dx*dx);
	printf("dt = %f, t = %f \n", dt, t);
	int timeSteps = (int) (t/dt + 1);
	int n = (int) (1.0/dx + 1);				// L = 1
	printf("%d\n", n);
	double alpha = dt/(dx*dx);
	printf("time steps = %d\n", timeSteps);

	double **u = createMatrix(3, n);
	double *exact = createVector(1, n); 		 // ------FIX THIS-------
	double **error = new double*[3];

	u[0][n-1] = 1;
	u[1][n-1] = 1;
	u[2][n-1] = 1;

	forwardEuler(u[0], alpha, timeSteps, n);
	backwardEuler(u[1], alpha, timeSteps, n);
	crank_nicolson(u[2], alpha, timeSteps, n);

	error[0] = relError(exact, u[0], n);
	error[1] = relError(exact, u[1], n);
	error[2] = relError(exact, u[2], n);


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

