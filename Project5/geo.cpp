#define _USE_MATH_DEFINES
#include "PDE.h"
#include "utils.h"	
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <cmath>

using namespace std;

void geo(double t, double h, string filename, int my_rank, int num_procs);

int main(int narg, char** argv) {
	int my_rank, num_procs;
	MPI_Init(&narg, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	double **upper_u, **middle_u, **lower_u;
	int **displs, **count;
	int upper_m, middle_m, lower_m;
	double x, y;
	double t = 1.0;
	double dt = 0.25*(h*h*h*h);
	int timeSteps = (int) (t/dt + 1);
	int n = (int) (1.0/h + 1);

	if (my_rank == 0) {
		count = createMatrix(2, num_procs);
		displs = createMatrix(2, num_procs);
		init(count[]);
	} else {
		MPI_Recv(&upper_m, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&y, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}



	MPI_Finalize();
	return 0;
}

void init(int *displs, int *count, int *my_m, int total_m) {

}

void geo(double t, double h, string filename, int my_rank, int num_procs){
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

		for (int i = 1; i < num_procs; i++) {
			displs[i] = displs[i-1] + count[i-1];
			if (n % num_procs == (num_procs - i)) {
				your_m++;
			} 
			if (i == num_procs-1) {
				your_m--;
				count[i] = your_m*n;
			} else {
				count[i] = (your_m-2)*n;
			}
			MPI_Send(&your_m, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&your_y, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			your_y += (your_m-2)*h;
		}
	} else {
		MPI_Recv(&m, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&y, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	double **u = createMatrix(m, n);
	double *exact = createVector(0.0, m*n);

	// edges filled with zero
	for (int i = 0; i < m; i++) {
		for (int j = 1; j < (n-1); j++) {
			x = j*h;
			u[i][j] = sin(pi*x)*sin(pi*y);							// ------FIX THIS-------
			exact[i*n + j] = u[i][j]*exp(-2*pi*pi*t);				// ------FIX THIS-------
		}
		y += h;
	}

	if (my_rank == num_procs-1) {
		for (int i = 0; i < n; i++) {
			exact[(m-1)*n + i] = 0.0;
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