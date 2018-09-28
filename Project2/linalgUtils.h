#include "armadillo"

using namespace arma;


mat copySymMatrixToArma(double **A, int n);

double** transpose(double **A, int n);

void sortEig(double *eigval, double **eigvec, int n);

