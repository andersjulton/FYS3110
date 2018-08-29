#include <fstream>
#include <string>

void printError(double *u, double *v, int n) {
	// write out difference between u-array and v-array

	double error;
	printf("\n");
	for (int i = 0; i < n; i++) {
		error = u[i] - v[i];
		printf("v = %.20f, u = %.20f, error = %.30f\n", v[i], u[i], error);
	}
	printf("\n");
}

void writeToFile(double *v, int n) {
	// write v-array to a txt-file

	std::ofstream myfile("n_"+std::to_string(n)+".txt");
	if (myfile.is_open()) {
		myfile << n << "\n";
		for (int i = 0; i < n; i++) {
			myfile << v[i] << "\n";
		}
	}
}

