
#include "armadillo"

using namespace arma;

vec armadillo_LU_solve(mat A, vec b);

mat copyMatrixToArma(double **A, int n);

int compareArmaMatrix(double **A, mat armaA, int n);

int compareArmaVector(double *b, vec armab, int n);
