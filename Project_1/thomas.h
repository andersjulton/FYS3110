#pragma once
#include "armadillo"

using namespace arma;

double *createVector(double value, int n);

double *solutionVector(int n);

void	generalTriDiaSolver(double *a, double *b, double *c, double *v, double *f, int n);

void	constTriDiaSolver(double a, double c, double *b, double *v, double *f, int n);

double *exactSolution(int n);

vec armadillo_LU_solve(mat A, vec b);