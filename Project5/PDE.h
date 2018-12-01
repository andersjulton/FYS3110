#pragma once

void forwardEuler(double *u, double alpha, int timeSteps, int n);

void backwardEuler(double *u, double alpha, int timeSteps, int n);

void crank_nicolson(double *u, double alpha, int timeSteps, int n);