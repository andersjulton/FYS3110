#include "PDE_2D.h"
#include "utils.h"	
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <cmath>

using namespace std;

void writeToFile(double t, double h, string filename, int my_rank, int num_procs);

int main(int narg, char** argv) {
	int my_rank, num_procs;
	MPI_Init(&narg, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);


	double t1 = 0.05;
	double t2 = 0.005;
	writeToFile(t1, 0.01, "t1_h_0.01", my_rank, num_procs);
	writeToFile(t2, 0.01, "t2_h_0.01", my_rank, num_procs);
	writeToFile(t1, 0.1, "t1_h_0.1", my_rank, num_procs);
	writeToFile(t2, 0.1, "t2_h_0.1", my_rank, num_procs);

	MPI_Finalize();
	return 0;
}

void writeToFile(double t, double h, string filename, int my_rank, int num_procs){
	double pi = 3.141592653589793238462643383279502884;
	double dt = 0.25*(h*h*h*h);					// h = dx = dy
	double alpha = dt/(h*h);
	int timeSteps = (int) (t/dt + 1);
	int n = (int) (1.0/h + 1);
	int m, mySize;
	double x = 0.0;
	double y;
	double *whole_u, *whole_error, *myPart, *error;
	int *count, *displs;

	if (my_rank == 0) {
		m = n/num_procs + 1;
		count = createVectorInt(0, num_procs);
		displs = createVectorInt(0, num_procs);
		count[0] = (m-2)*n;
		y = 0.0;

		double your_y = (m-2)*h;
		int your_m = m + 1;

		printf("count = %d\n", count[0]);
		for (int i = 1; i < num_procs; i++) {
			displs[i] = displs[i-1] + count[i-1];
			printf("disp = %d\n", displs[i]);
			if (n % num_procs == (num_procs - i)) {
				your_m++;
			} 
			if (i == num_procs-1) {
				your_m--;
				count[i] = your_m*n;
			} else {
				count[i] = (your_m-2)*n;
			}
			printf("count = %d\n", count[i]);
			MPI_Send(&your_m, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&your_y, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			your_y += (your_m-2)*h;
		}
		printf("h = %f\n", h);
		printf("n = %d, n x n = %d, y = %f\n", n, n*n, your_y+h);
	} else {
		MPI_Recv(&m, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&y, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}


	printf("my rank is %d, m = %d, y = %f\n", my_rank, m, y);
	double **u = createMatrix(m, n);
	double *exact = createVector(0.0, m*n);

	// edges filled with zero
	for (int i = 0; i < m; i++) {
		for (int j = 1; j < (n-1); j++) {
			x = j*h;
			// using-ish eq. from page 314, just needed something to test, no idea if it's correct
			// using L = 1 & n = 1
			u[i][j] = sin(pi*x)*sin(pi*y);							// ------FIX THIS-------
			exact[i*n + j] = u[i][j]*exp(-2*pi*pi*t);				// ------FIX THIS-------
			//u[i][j] = my_rank + 1;
		}
		y += h;
	}

	if (my_rank == num_procs-1) {
		for (int i = 0; i < n; i++) {
			exact[(m-1)*n + i] = 0.0;
			u[(m-1)][i] = 0.0;
		}
	}
	forwardEuler(u, alpha, timeSteps, m, n, my_rank, num_procs-1);


	if (my_rank != num_procs-1) {
		mySize = (m-2)*n;
		myPart = copyMatrixTo1D(u, m-2, n);
		error = absError(exact, myPart, mySize);
	} else {
		mySize = m*n;
		myPart = copyMatrixTo1D(u, m, n);
		error = absError(exact, myPart, mySize);
	}
	if (my_rank == 0) {
		whole_u = createVector(0.0, n*n);
		whole_error = createVector(0.0, n*n);
	}
	printf("rank (%d) my size is %d\n", my_rank, mySize);

	MPI_Gatherv(myPart, mySize, MPI_DOUBLE, whole_u, count, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(error, mySize, MPI_DOUBLE, whole_error, count, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (my_rank == 0) {
		doubleMatrixToFile(whole_u , n, n, filename);
		doubleMatrixToFile(whole_error , n, n, filename + "_error");
		delete[] whole_u;
		delete[] whole_error;
		delete[] count;
		delete[] displs;
	}
	
	deleteMatrix(u, m);
	delete[] exact;
	delete[] error;
	delete[] myPart;
}