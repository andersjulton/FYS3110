#include <fstream>
#include "func.h"
#include "utils.h"
#include "armadillo"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <string>

using namespace std;
using namespace arma;

// Solving Ax = b
vec armadillo_LU_solve(mat A, vec b) {
	mat L, U;
	lu(L, U, A);
	// only for optimization 
	mat L_prime = trimatl(L);
	mat U_prime = trimatu(U);
	//To solve LUx = b, you first solve Ly = b and then Ux = y.
	vec y, x;
	y = solve(L_prime, b);
	x = solve(U_prime, y);
	return x;
}

int testOrthogonality(double **A){


	return 0;
}