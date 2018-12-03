#pragma once

void forwardEuler(double **u, double alpha, int timeSteps, int m, int n);

void forwardEuler(double **u, double alpha, int timeSteps, int m, int n, int my_rank, int final_rank);