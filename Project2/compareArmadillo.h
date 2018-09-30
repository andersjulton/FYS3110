#pragma once
#include "armadillo"

using namespace arma;

mat copySymMatrixToArma(double **A, int n);

int compareArmaMatrix(double **A, mat armaA, int n);

int compareArmaVector(double *b, vec armab, int n);


