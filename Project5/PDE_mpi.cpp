#include "triDiaSolver.h"
#include "utils.h"
#include <mpi.h>

using namespace std;

// parallel implementation of the explicit forward Euler algorithm in 2D
void forwardEuler(double **u, double alpha, int timeSteps, int m, int n, int my_rank, int final_rank) {
	double **u_new = copyMatrix(u, m, n);
	double **temp;

	for (int time = 0; time < timeSteps; time++) {	
		for (int i = 1; i < (m-1); i++) {
			for (int j = 1; j < (n-1); j++) {
				u_new[i][j] = u[i][j] + alpha*(u[i+1][j] + u[i-1][j]+ u[i][j+1] + u[i][j-1] - 4*u[i][j]);
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
