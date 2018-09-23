#include <string>
#include <cmath>
#include <fstream>

using namespace std; 	// Oh, no!

//  Allocating space for a vector and fill elements with a constant value
double* createVector(double value, int n) {
	double *vector;
	vector = new double[n];
	for (int i = 0; i < n; i++) {
		vector[i] = value;
	}
	return vector;
}

// returns the largest relative error comparing two arrays
double maxError(double *expected, double *computed, int n) {
	double max, max_i;
	max = fabs((expected[0] - computed[0])/expected[0]);
	for (int i = 1; i < n; i++) {
		max_i = fabs((expected[i] - computed[i])/expected[i]);
		if (max < max_i) {
			max = max_i;
		}
	}
	return max;
}

double maxEpsilon(double *expected, double *computed, int n) {
	double max = maxError(expected, computed, n);
	return log10(max);
}

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

// write single array to a txt-file
void arrayToFile(double *v , int n, string filename) {
	ofstream myfile(filename + ".txt");
	if (myfile.is_open()) {
		myfile << n << "\n";
		for (int i = 0; i < n; i++) {
			myfile << v[i] << "\n";
		}
	}
}