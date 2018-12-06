#include "triDiaSolver.h"
#include "utils.h"
#include <mpi.h>

using namespace std;

void forwardEuler(double **u, double alpha, int timeSteps, int m, int n, int my_rank, int final_rank) {
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

		// sending new edges
		if ((my_rank % 2) == 0) {
			if (my_rank != 0) {
				MPI_Send(u[1], n, MPI_DOUBLE, (my_rank-1), 0, MPI_COMM_WORLD);
				MPI_Recv(u[0], n, MPI_DOUBLE, (my_rank-1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			if (my_rank != final_rank) {
				MPI_Recv(u[m-1], n, MPI_DOUBLE, (my_rank+1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Send(u[m-2], n, MPI_DOUBLE, (my_rank+1), 0, MPI_COMM_WORLD);
			}
		} else {
			if (my_rank != final_rank) {
				MPI_Recv(u[m-1], n, MPI_DOUBLE, (my_rank+1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Send(u[m-2], n, MPI_DOUBLE, (my_rank+1), 0, MPI_COMM_WORLD);
			}
			// rank != 0 -> true
				MPI_Send(u[1], n, MPI_DOUBLE, (my_rank-1), 0, MPI_COMM_WORLD);
				MPI_Recv(u[0], n, MPI_DOUBLE, (my_rank-1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}	
	}
	deleteMatrix(u_new, m);
}
