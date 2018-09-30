#pragma once

int jacobi(double **A, double **R, int n, double epsilon);

int* getMaxInRow(double **A, int n);

double maxOffDiag(double **A, int *indexOfMax, int *k, int *l, int n);

void updateMaxInRow(double **A, int *indexOfMax, int k, int l, int n);

void rotate(double **A, double **R, int k, int l, int n);

double *bisect(double *off, double *dia, int n, double epsilon);

int isolate(double xmin, double xmax, double *off, double *dia, double *eig, int count, int n, double epsilon);

double extract(double xmin, double xmax, double *off, double *dia, int n, double epsilon);

double *getInterval(double *off, double *dia, int n);

int getSignChange(double *a, double *b, int n, double lambda);

