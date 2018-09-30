#pragma once

double** transpose(double **A, int n);

void sortEig(double *eigval, double **eigvec, int n);

int testEig(double *eigval, double **eigvec, double **A, int n);

int testOrthogonality(double **A, int n);

double **transpose(double **A, int n);

double analyticConvergenceRate(int n, double eps, double sumOff);

void normalize(double *v, double *u, int n);

void extractEigenVec(double **A, double *u, int index, int n);