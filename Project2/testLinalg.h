#include "armadillo"

using namespace arma;


int compareArmaMatrix(double **A, mat armaA, int n);

int compareArmaVector(double *b, vec armab, int n);

int testVector(double *u, double *v, int n);

int testEig(double *eigval, double **eigvec, double **A, int n);

int testOrthogonality(double **A, int n);
