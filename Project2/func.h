#pragma once

double *solutionVector(int n);

void	generalTriDiaSolver(double *a, double *b, double *c, double *v, double *f, int n);

void	constTriDiaSolver(double a, double b, double c, double *v, double *f, int n);

double *exactSolution(int n);