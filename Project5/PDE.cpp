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

	double lowerBound = u[0];
	double upperBound = u[n-1];

	double *r = createVector(0, n);
	r[0] = lowerBound;
	r[n-1] = upperBound;

	for (int i = 0; i < timeSteps; i++) {		// time
		// Crank Nicolson method for right hand side
		for (int j = 1; j < (n-1); j++) {
			r[j] = 2*u[j]*(1 - alpha) + alpha*(u[j-1] + u[j+1]);
		}
		triDiaSolver(offDia, dia, offDia, u, r, n);
	}
	delete[] r;
}