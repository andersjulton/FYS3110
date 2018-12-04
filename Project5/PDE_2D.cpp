#include "triDiaSolver.h"
#include "utils.h"
#include <mpi.h>

using namespace std;

//u[y][x]
void forwardEuler(double **u, double alpha, int timeSteps, int m, int n) {
	double **u_new = copyMatrix(u, m, n);
	double **temp;


	for (int i = 0; i < timeSteps; i++) {		// time
		for (int j = 1; j < (m-1); j++) {
			for (int k = 1; k < (n-1); k++) {
				u_new[j][k] = u[j][k] + alpha*(u[j+1][k] + u[j-1][k]+ u[j][k+1] + u[j][k-1] - 4*u[j][k]);
			}
		}
		temp = u;
		u = u_new;
		u_new = temp;
		
	}
	deleteMatrix(u_new, m);
}



