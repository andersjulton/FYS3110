#pragma once

// 1D:
void forwardEuler(double *u, double alpha, int timeSteps, int n);

void backwardEuler(double *u, double alpha, int timeSteps, int n);

void crank_nicolson(double *u, double alpha, int timeSteps, int n);


// 2D:
void forwardEuler(double **u, double alpha, int timeSteps, int m, int n); 

void jacobi(double **u, double alpha, int timeSteps, int m, int n);

// 2D & parallel:
void forwardEuler(double **u, double alpha, int timeSteps, int m, int n, int my_rank, int final_rank);

void jacobi(double **u, double alpha, int timeSteps, int m, int n, int my_rank, int final_rank);