#pragma once

int jacobi(double **A, int n);

void rotate(double **A, double **R, int k, int l, int n);

double maxOffDiag (double **A, int *k, int *l, int n);

void bisect(double *a, double *b, int n);

int getSignChange(double *a, double *b, int n, double lambda);

double *getInterval(double *a, double *b, int n);