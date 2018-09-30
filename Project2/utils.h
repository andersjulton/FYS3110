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

void intArrayToFile(int *v , int n, std::string filename, bool zeroPadding = false);

void doubleArrayToFile(double *v , int n, std::string filename, bool zeroPadding = false);
