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

//  Allocating space for a m x n matrix and fill elements with 0
double **createMatrix(int m, int n) {
	double **mat;
	mat = new double*[m];
	for(int i = 0; i < m; i++) {
		mat[i] = new double[n];
		for(int j = 0; j < n; j++) {
			mat[i][j] = 0.0;
		}
	}
	return mat;
}

// delete 2D-array
void deleteMatrix(double **mat, int n) {
	for (int i = 0; i < n; i++) {
		delete[] mat[i];
	}
	delete[] mat;
}

// Allocating space for a n x n diagonal matrix and fill the diagonal with d
double **createDiaMatrix(double d, int n) {
	// two seperate loops to avoid if-tests
	double **mat;
	mat = createMatrix(n, n);
	for(int i = 0; i < n; i++) {
		mat[i][i] = d;
	}
	return mat;
}

// Create a tridiagonal matrix
double **createTriDiaMatrix(double off_value, double d_value, int n) {
	// two seperate loops to avoid if-tests
	double **mat;
	mat = createMatrix(n, n);
	mat[0][0] = d_value;
	mat[0][1] = off_value;
	for (int i = 1; i < (n - 1); i++) {
		mat[i][i-1] = off_value;
		mat[i][i] = d_value;
		mat[i][i+1] = off_value;
	}
	mat[n-1][n-2] = off_value;
	mat[n-1][n-1] = d_value;
	return mat;
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