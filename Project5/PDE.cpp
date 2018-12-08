#include "triDiaSolver.h"
#include "utils.h"

using namespace std;

// implementation of the explicit forward Euler algorithm in 1D
void forwardEuler(double *u, double alpha, int timeSteps, int n) {
	double u_prev, u_i;		// avoiding several vectors

	for (int time = 0; time < timeSteps; time++) {
		u_prev = u[0];
		for (int i = 1; i < n-1; i++) {
			u_i = u[i];
			u[i] = u_i + alpha*(u_prev + u[i+1] - 2*u_i);
			u_prev = u_i;
		}
	}
}

// implementation of the implicit backward Euler algorithm in 1D
void backwardEuler(double *u, double alpha, int timeSteps, int n) {
	double offDia = -alpha;
	double dia = 1 + 2*alpha;

	double *u_new = createVector(0, n);

	for (int time = 0; time < timeSteps/2; time++) {
		// avoiding copying between vectors by doing two time-iterations at once
		triDiaSolver(offDia, dia, offDia, u_new, u, n);
		triDiaSolver(offDia, dia, offDia, u, u_new, n);
	}
	if (timeSteps % 2 == 1) {
		triDiaSolver(offDia, dia, offDia, u_new, u, n);
		for (int i = 1; i < (n-1); i++) {
			u[i] = u_new[i];
		}
	}
	delete[] u_new;
}

// implementation of the implicit Crank-Nicolson scheme in 1D
void crank_nicolson(double *u, double alpha, int timeSteps, int n) {
	double offDia = -alpha;
	double dia = 2 + 2*alpha;

	double *r = createVector(0, n);
	r[0] = u[0];
	r[n-1] = u[n-1];

	for (int time = 0; time < timeSteps; time++) {	
		// Crank Nicolson method for right hand side
		for (int j = 1; j < (n-1); j++) {
			r[j] = 2*u[j]*(1 - alpha) + alpha*(u[j-1] + u[j+1]);
		}
		triDiaSolver(offDia, dia, offDia, u, r, n);
	}
	delete[] r;
}

// implementation of the explicit forward Euler algorithm in 2D
void forwardEuler(double **u, double alpha, int timeSteps, int m, int n) {
	double **u_new = copyMatrix(u, m, n);
	double **temp;

	for (int time = 0; time < timeSteps; time++) {	
		for (int i = 1; i < (m-1); i++) {
			for (int j = 1; j < (n-1); j++) {
				u_new[i][j] = u[i][j] + alpha*(u[i+1][j] + u[i-1][j]+ u[i][j+1] + u[i][j-1] - 4*u[i][j]);
			}
		}
		temp = u;
		u = u_new;
		u_new = temp;
	}
	deleteMatrix(u_new, m);
}