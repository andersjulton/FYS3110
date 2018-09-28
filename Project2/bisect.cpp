#include <cmath>
#include <stdlib.h>
#include "utils.h"
#include "bisect.h"
#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>

using namespace std;

double *Gershgorin(double *off, double *dia, int n) {
	double *interval = getInterval(off, dia, n);
	double xmin = interval[0];
	double xmax = interval[1];
	printf("X min = %.5f, max = %.5f\n", xmin, xmax);
	int signMin = getSignChange(off, dia, n, xmin);
	int signMax = getSignChange(off, dia, n, xmax);
	int numOfEig = signMax - signMin;
	printf("S min = %d, max = %d\n", signMin, signMax);
	printf("numEig = %d\n", numOfEig);
	double *eig = createVector(0, n);
	isolate(xmin, xmax, off, dia, eig, 0, n);
	return eig;
}

int isolate(double xmin, double xmax, double *off, double *dia, double *eig, int count, int n) {
	int signMin = getSignChange(off, dia, n, xmin);
	int signMax = getSignChange(off, dia, n, xmax);
	int numOfEig = signMax - signMin;
	//printf("X min = %.5f, max = %.5f\n", xmin, xmax);
	//printf("S min = %d, max = %d\n", signMin, signMax);
	if (numOfEig == 0) {
		return count;
	} else if (numOfEig == 1) {
		//printf("X min = %.5f, max = %.5f\n", xmin, xmax);
		//printf("S min = %d, max = %d\n", signMin, signMax);
		eig[count] = bisect(xmin, xmax, off, dia, n);
		count++;
		return count;
	} else {
		if (numOfEig < 0) {
			printf("HALLA\n");
			return count;
		}
		double xmiddle = xmin + (xmax - xmin)/2.0;
		count = isolate(xmin, xmiddle, off, dia, eig, count, n);
		count = isolate(xmiddle, xmax, off, dia, eig, count, n);
	}
	return count;
}


// call when final interval found
double bisect(double xmin, double xmax, double *off, double *dia, int n) {
	double eps = 1e-13;
	int maxIter = 200;
	int iter = 0;
	int signMin, signMax, numOfEig;
	double x;
	while (fabs(xmax - xmin) > eps && iter < maxIter) {
		x = xmin + (xmax - xmin)/2.0;
		signMin = getSignChange(off, dia, n, xmin);
		signMax = getSignChange(off, dia, n, x);
		numOfEig = signMax - signMin;
		if (numOfEig == 1) {
			xmax = x;
		} else {
			xmin = x;
		}
		iter++;
	}
	return xmin + (xmax - xmin)/2.0;
}

// Find max value of interval by Gershgorin’s Theorem
double *getInterval(double *off, double *dia, int n) {
	double xmax, xmin, h;
	double *interval;
	interval = createVector(0, 2);
	xmin = dia[0] - fabs(off[0]);
	xmax = dia[0] + fabs(off[0]);
	for (int i = 1; i < (n-1); i++) {
		h = fabs(off[i-1]) + fabs(off[i]);
		if ((dia[i] - h) < xmin) {
			xmin = dia[i] - h;
		}
		if ((dia[i] + h) > xmax) {
			xmax = dia[i] + h;
		}
	}
	h = fabs(off[n-2]);
	if ((dia[n-1] - h) < xmin) {
		xmin = dia[n-1] - h;
	}
	if ((dia[n-1] + h) > xmax) {
		xmin = dia[n-1] + h;
	}
	interval[0] = xmin;
	interval[1] = xmax;
	return interval;
}


int getSignChange(double *off, double *dia, int n, double lambda) {
	int count = 0;
	double eps = 1e-12;
	double q;
	double pre = dia[0] - lambda;
	if (pre < 0) {
		count++;
	} else if (pre == 0.0) {
		pre = eps;
	}
	for (int i = 1; i < n; i++) {
		q = (dia[i] - lambda) - (off[i-1]*off[i-1])/pre;
		if (q < 0) {
			count++;
		} else if (q == 0.0) {
			q = eps;
		}
		pre = q;
	}
	return count;
}

