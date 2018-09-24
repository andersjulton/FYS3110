
#include "armadillo"

using namespace arma;

vec armadillo_LU_solve(mat A, vec b);

double maxErrorDiaSolver(int n);

void compareTime(int a_value, int b_value, int c_value, int n);

void compareTimeArmadillo(int n);
