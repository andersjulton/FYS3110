#include "jacobi_bisect.h"
#include "utils.h"
#include "linalgUtils.h"
#include "armadillo"
#include "compareArmadillo.h"
#include <string>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <chrono>

//using namespace arma;
using namespace std;

void oneElectron(int n, double epsilon, double rho_N);
void twoElectrons(int n, double omega, std::string filename, double epsilon, bool writeToFile);
void bucklingBeam(int n, double epsilon);
void takeTime(int n);

int main() {
	
	bucklingBeam(5, 1e-8);
	bucklingBeam(10, 1e-8);
	bucklingBeam(50, 1e-8);
	bucklingBeam(200, 1e-8);

	oneElectron(10, 1e-5, 5);
	oneElectron(100, 1e-5, 5);
	oneElectron(200, 1e-5, 5);
	oneElectron(400, 1e-5, 5);

	int n = 400;
	std::string filename1 = "wave_1";
	std::string filename2 = "wave_2";
	std::string filename3 = "wave_3";	
	std::string filename4 = "wave_4";

	bool print = true;
	twoElectrons(n, 0.01, filename1, 1e-5, print);
	twoElectrons(n, 0.5, filename2, 1e-5, print);
	twoElectrons(n, 1.0, filename3, 1e-5, print);
	twoElectrons(n, 5.0, filename4, 1e-5, print);
	//system("PAUSE");

	printf("\n");
	takeTime(10);
	takeTime(50);
	takeTime(100);
	takeTime(200);
	takeTime(300);
	takeTime(400);
	takeTime(500);

	return 0;
}

void bucklingBeam(int n, double epsilon) {
	printf("\nCalculating eigenvalues for the buckling beam for n = %d \n\n", n);
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

	double *eig = diagToVector(A, n);
	double **E = transpose(R, n);
	sortEig(eig, E, n);

	double *off = createVector(a, n);
	double *dia = createVector(d, n);
	double *l = bisect(off, dia, n, epsilon);

	double lambda;
	for (int i = 0; i < 4; i++) {
		lambda = d + 2*a*cos((i+1)*pi/(n + 1));
		printf("Jacobi = %5.4f  Bisect = %5.4f  Analytical = %5.4f\n", eig[i], l[i], lambda);
	}
	deleteMatrix(E, n);
	deleteMatrix(A, n);
	deleteMatrix(R, n);
	delete[] eig;
	delete[] u;
	delete[] off;
	delete[] dia;
	delete[] l;
}

void oneElectron(int n, double epsilon, double rho_N) {
	printf("\nCalculating eigenvalues for one electron for n = %d and rho_max = %2.1f\n\n", n, rho_N);
	int iterations;
	double rho_0, rho_i, h;
	double **A, **R;
	double *u;
	u = createVector(0, n);
	R = createDiaMatrix(1, n);
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

	double *eig = diagToVector(A, n);
	double **E = transpose(R, n);
	double lambda;
	sortEig(eig, E, n);
	
	for (int i = 0; i < 4; i++) {
		lambda = 3.0 + 4.0*i;
		printf("Calculated = %2.5f. Analytical = %1.0f  Error = %5.5f\n", eig[i], lambda, fabs((eig[i] - lambda)/lambda));
	}
	
	deleteMatrix(A, n);
	deleteMatrix(R, n);
	deleteMatrix(E, n);
	delete[] eig;
	delete[] u;
}

void twoElectrons(int n, double omega, std::string filename, double epsilon, bool writeToFile) {
	printf("\nCalculating eigenvalues for two electrons for n = %d and omega = %3.2f\n", n, omega);
	int iterations;
	double rho_0, rho_N, rho_i, h;
	double **A, **R;
	double *u, *un;
	u = createVector(0, n);
	un = createVector(0, n);
	R = createDiaMatrix(1, n);
	rho_N = 5;
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

	printf("----Eigenvalue in ground state = %5.5f\n", eig[0]);

	extractEigenVec(E, u, 0, n);
	normalize(u, un, n);

	if (writeToFile) {
		doubleArrayToFile(u, n, filename, true);
		doubleArrayToFile(un, n, filename+"n", true);		
	}
	
	deleteMatrix(A, n);
	deleteMatrix(R, n);
	deleteMatrix(E, n);
	delete[] eig;
	delete[] u;
	delete[] un;
}

void takeTime(int n) {
	chrono::duration<double> elapsed;

	double *timeJacobi = new double[5];
	double *timeBisect = new double[5];
	double *timeArmadillo = new double[5];

	double offDiaValue = -1.0;
	double diaValue = 2.0;

	double **A, **R;
	mat armaA;
	vec eig;
	double *off = createVector(offDiaValue, n);
	double *dia = createVector(diaValue, n);
	printf("n = %d\n", n);
	// timing 5 times
	for (int i = 0; i < 5; i++) {
		A = createTriDiaMatrix(offDiaValue, diaValue, n);
		armaA = copySymMatrixToArma(A, n);
		R = createDiaMatrix(1, n);
		// time Jacobi
		auto begin = chrono::high_resolution_clock::now();
		jacobi(A, R, n, 1e-8);
		auto end = chrono::high_resolution_clock::now();
		elapsed = (end - begin);
		timeJacobi[i] = (double)elapsed.count();
		// time Bisect
		begin = chrono::high_resolution_clock::now();
		bisect(off, dia, n, 1e-8);
		end = chrono::high_resolution_clock::now();
		elapsed = (end - begin);
		timeBisect[i] = (double)elapsed.count();
		// time armadillo
		begin = chrono::high_resolution_clock::now();
		eig = eig_sym(armaA);
		end = chrono::high_resolution_clock::now();
		elapsed = (end - begin);
		timeArmadillo[i] = (double)elapsed.count();
		deleteMatrix(A, n);
		deleteMatrix(R, n);
	}
	// Finding the median
	sort(timeJacobi, timeJacobi + 5);
	printf("Time it took for Jacobi: %g s\n", timeJacobi[2]);

	sort(timeBisect, timeBisect + 5);
	printf("Time it took for Bisect: %g s\n", timeBisect[2]);

	sort(timeArmadillo, timeArmadillo + 5);
	printf("Time it took for armadillo: %g s\n", timeArmadillo[2]);

	printf("Ratio (armadillo/jacobi): %g \n", timeArmadillo[2]/timeJacobi[2]);
	printf("Ratio (armadillo/bisect): %g \n", timeArmadillo[2]/timeBisect[2]);
	printf("Ratio (bisect/jacobi): %g \n", timeBisect[2]/timeJacobi[2]);
	printf("Ratio (jacobi/armadillo): %g \n", timeJacobi[2]/timeArmadillo[2]);
	printf("Ratio (bisect/armadillo: %g \n", timeBisect[2]/timeArmadillo[2]);
	printf("Ratio (jacobi/bisect): %g \n", timeJacobi[2]/timeBisect[2]);

	delete[] off;
	delete[] dia;
	delete[] timeBisect;
	delete[] timeJacobi;
	delete[] timeArmadillo; 
}