#include "triDiaSolver.h"
#include "utils.h"

using namespace std;

void forwardEuler(double *u, double alpha, int timeSteps, int n) {
	double u_prev, u_j;		// avoiding several vectors

	for (int i = 0; i < timeSteps; i++) {		// time
		u_prev = u[0];
		for (int j = 1; j < n-1; j++) {
			u_j = u[j];
			u[j] = u_j + alpha*(u_prev + u[j+1] - 2*u_j);
			u_prev = u_j;
		}
	}
}

void backwardEuler(double *u, double alpha, int timeSteps, int n) {
	double offDia = -alpha;
	double dia = 1 + 2*alpha;

	double *temp = createVector(0, n);

	for (int i = 0; i < timeSteps/2; i++) {		// time
		// avoiding copying between vectors by doing two time-iterations at once
		triDiaSolver(offDia, dia, offDia, temp, u, n);
		triDiaSolver(offDia, dia, offDia, u, temp, n);
	}
	if (timeSteps % 2 == 1) {
		triDiaSolver(offDia, dia, offDia, temp, u, n);
		for (int i = 1; i < (n-1); i++) {
			u[i] = temp[i];
		}
	}
	delete[] temp;
}

void crank_nicolson(double *u, double alpha, int timeSteps, int n) {
	double offDia = -alpha;
	double dia = 2 + 2*alpha;

	double *r = createVector(0, n);
	r[0] = u[0];
	r[n-1] = u[n-1];

	for (int i = 0; i < timeSteps; i++) {		// time
		// Crank Nicolson method for right hand side
		for (int j = 1; j < (n-1); j++) {
			r[j] = 2*u[j]*(1 - alpha) + alpha*(u[j-1] + u[j+1]);
		}
		triDiaSolver(offDia, dia, offDia, u, r, n);
	}
	delete[] r;
}

// 2D
void forwardEuler(double **u, double alpha, int timeSteps, int m, int n) {
	double **u_new = copyMatrix(u, m, n);
	double **temp;

	for (int i = 0; i < timeSteps; i++) {		// time
		for (int j = 1; j < (m-1); j++) {
			for (int k = 1; k < (n-1); k++) {
				u_new[j][k] = u[j][k] + alpha*(u[j+1][k] + u[j-1][k]+ u[j][k+1] + u[j][k-1] - 4*u[j][k]);
			}
		}
		temp = u;
		u = u_new;
		u_new = temp;
	}
	deleteMatrix(u_new, m);
}