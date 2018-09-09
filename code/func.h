#pragma once
#include "armadillo"

using namespace arma;

void	printError(double *u, double *v, int n);

void	writeToFile(double *vec, int n);

double *createVector(double value, int n);

double *solutionVector(int n);

void	generalTriDiaSolver(double *a, double *b, double *c, double *v, double *f, int n);

void	constTriDiaSolver(double a, double c, double *b, double *v, double *f, int n);

double *exactSolution(int n);

vec armadillo_LU_solve(mat A, vec b);

double maxError(double *u, double *v, int n);

double maxErrorDiaSolver(int n);

void compareTime(int a_value, int b_value, int c_value, int n);

void compareTimeArmadillo(int n);