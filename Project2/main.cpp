#include "jacobi.h"
#include "utils.h"
#include <string>
#include <algorithm>
#include <iostream>
#include <cmath>

void oneElectron(int n, double epsilon);
void twoElectrons(int n, double omega, std::string filename, double epsilon);
void bucklingBeam(int n, double epsilon);

int main() {
	for (int i = 5; i < 20; i++) {
		bucklingBeam(i, 1e-5);
	}
	
	
	system("PAUSE");
	/*
	int n = 200;
	double omega = 0.5;
	std::string filename1 = "wave_21";
	std::string filename2 = "wave_22";
	std::string filename3 = "wave_23";	
	std::string filename4 = "wave_24";
	//oneElectron(n);
	twoElectrons(n, 0.01, filename1);
	twoElectrons(n, 0.5, filename2);
	twoElectrons(n, 1.0, filename3);
	twoElectrons(n, 5.0, filename4);*/
	return 0;
}

void bucklingBeam(int n, double epsilon) {
	int iterations;
	double **A, **R;
	double *u;
	double h, d, a, maxInt;
	double pi = 3.14159265359;
	u = createVector(0, n);
	R = createDiaMatrix(1, n);
	h = 1.0/n;
	a = -1.0/(h*h);
	maxInt = analyticConvergenceRate(n, epsilon, 2*(n-1)*a*a);
	
	d = 2.0/(h*h);
	A = createTriDiaMatrix(a, d, n);
	iterations = jacobi(A, R, n, epsilon);
	printf("\nIterations \n");
	printf("n = %d  Actual iterations = %d  Analytical iterations = %5.0f\n", n, iterations, maxInt);
	double *eig = diagToVector(A, n);
	double **E = transpose(R, n);
	double lambda;
	sortEig(eig, E, n);
	printf("\nEigenvalues \n");
	for (int i = 0; i < n; i++) {
		lambda = d + 2*a*cos((i+1)*pi/(n + 1));
		
		//printf("Calculated = %5.4f     Actual = %5.4f\n", eig[i], lambda);
	}
	deleteMatrix(E, n);
	deleteMatrix(A, n);
	deleteMatrix(R, n);
	delete[] eig;
	delete[] u;
}

void oneElectron(int n, double epsilon) {
	int iterations;
	double rho_0, rho_N, rho_i, h;
	double **A, **R;
	double *u;
	u = createVector(0, n);
	R = createDiaMatrix(1, n);
	rho_N = 10.0;
	rho_0 = 0;
	h = (rho_N - rho_0)/n;
	A = createMatrix(n, n);
	A[0][0] = 2.0/(h*h);
	A[0][1] = -1.0/(h*h);
	for (int i = 1; i < (n - 1); i++) {
		rho_i = (i+1)*h;
		A[i][i-1] = -1.0/(h*h);
		A[i][i] = 2.0/(h*h) + rho_i*rho_i;
		A[i][i+1] = -1.0/(h*h);
	}
	rho_i = (n - 1);
	A[n-1][n-2] = -1.0/(h*h);
	A[n-1][n-1] = 2.0/(h*h) + rho_i*rho_i;

	iterations = jacobi(A, R, n, epsilon);
	std::cout << iterations << '\n';

	double *eig = diagToVector(A, n);
	double **E = transpose(R, n);
	sortEig(eig, E, n);

	for (int j = 0; j < n; j++) {
		u[j] = eig[0]*E[0][j];
	}
	//arrayToFile(u, n, "wave_1");
	deleteMatrix(A, n);
	deleteMatrix(R, n);
	deleteMatrix(E, n);
	delete[] eig;
	delete[] u;
}

void twoElectrons(int n, double omega, std::string filename, double epsilon) {
	int iterations;
	double rho_0, rho_N, rho_i, h;
	double **A, **R;
	double *u;
	u = createVector(0, n);
	R = createDiaMatrix(1, n);
	rho_N = 10.0;
	rho_0 = 0;
	h = (rho_N - rho_0)/n;
	A = createMatrix(n, n);
	A[0][0] = 2.0/(h*h);
	A[0][1] = -1.0/(h*h);
	for (int i = 1; i < (n - 1); i++) {
		rho_i = (i+1)*h;
		A[i][i-1] = -1.0/(h*h);
		A[i][i] = 2.0/(h*h) + rho_i*rho_i*omega*omega + 1.0/rho_i;
		A[i][i+1] = -1.0/(h*h);
	}
	rho_i = (n - 1);
	A[n-1][n-2] = -1.0/(h*h);
	A[n-1][n-1] = 2.0/(h*h) + rho_i*rho_i*omega*omega + 1.0/rho_i;

	iterations = jacobi(A, R, n, epsilon);
	double *eig = diagToVector(A, n);
	double **E = transpose(R, n);
	sortEig(eig, E, n);


	for (int j = 0; j < n; j++) {
		u[j] = eig[0]*E[0][j];
	}
	arrayToFile(u, n, filename);
	deleteMatrix(A, n);
	deleteMatrix(R, n);
	deleteMatrix(E, n);
	delete[] eig;
	delete[] u;
}