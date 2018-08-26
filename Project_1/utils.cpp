#include <fstream>
#include "utils.h"

// write out difference between u-array and v-array
void printError(double *u, double *v, int n) {
	double error;
	printf("\n");
	for (int i = 0; i < n; i++) {
		error = u[i] - v[i];
		printf("v = %.20f, u = %.20f, error = %.30f\n", v[i], u[i], error);
	}
	printf("\n");
}

// write v-array to a txt-file
void writeToFile(double *v, int n) {
	std::ofstream myfile("n_"+std::to_string(n)+".txt");
	if (myfile.is_open()) {
		myfile << n << "\n";
		for (int i = 0; i < n; i++) {
			myfile << v[i] << "\n";
		}
	}
}

