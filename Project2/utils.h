#pragma once
#include <string>


double* diagToVector(double **A, int n);

double* createVector(double value, int n);

double **createMatrix(int m, int n);

double **createDiaMatrix(double d, int n);

double **createTriDiaMatrix(double off_value, double d_value, int n);

void deleteMatrix(double **mat, int n);

double maxError(double *expected, double *computed, int n);

double maxEpsilon(double *expected, double *computed, int n);

void printError(double *u, double *v, int n);

void arrayToFile(double *v , int n, std::string filename, bool zeroPadding = false);

void sortEig(double *eigval, double **eigvec, int n);

double **transpose(double **A, int n);

double analyticConvergenceRate(int n, double eps, double sumOff);

void normalize(double *v, double *u, int n);

void extractEigenVec(double **A, double *u, int index, int n);