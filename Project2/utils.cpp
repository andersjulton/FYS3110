#include <string>
#include <cmath>
#include <fstream>


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

// Allocating space for a n x n diagonal matrix and fill the diagonal with d
double **createDiaMatrix(double d, int n) {
	double **mat;
	mat = createMatrix(n, n);
	// two seperate loops to avoid if-tests
	for(int i = 0; i < n; i++) {
		mat[i][i] = d;
	}
	return mat;
}

// Create a tridiagonal matrix
double **createTriDiaMatrix(double off_value, double d_value, int n) {
	double **mat;
	mat = createMatrix(n, n);
	mat[0][0] = d_value;
	mat[0][1] = off_value;
	// two seperate loops to avoid if-tests
	for (int i = 1; i < (n - 1); i++) {
		mat[i][i-1] = off_value;
		mat[i][i] = d_value;
		mat[i][i+1] = off_value;
	}
	mat[n-1][n-2] = off_value;
	mat[n-1][n-1] = d_value;
	return mat;
}

// delete 2D-array
void deleteMatrix(double **mat, int n) {
	for (int i = 0; i < n; i++) {
		delete[] mat[i];
	}
	delete[] mat;
}

// Represent diagonal matrix A as a vector
double* diagToVector(double **A, int n) {
	double *a = new double[n];
	for (int i = 0; i < n; i++) {
		a[i] = A[i][i];
	}
	return a;
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

// returns log10 of the max relative error
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
void arrayToFile(double *v , int n, std::string filename, bool zeroPadding = false) {
	std::ofstream myfile(filename + ".txt");
	if (myfile.is_open()) {
		if (zeroPadding) {
			myfile << n+2 << "\n";
			myfile << 0.0 << "\n";
		}
		else {
			myfile << n << "\n";
		}
		for (int i = 0; i < n; i++) {
			myfile << v[i] << "\n";
		}
		if (zeroPadding) {
			myfile << 0.0 << "\n";
		}
	}
}

void sortEig(double *eigval, double **eigvec, int n) {
	int j;
	double current, past;
	double *currentVec, *pastVec;
	for (int i = 1; i < n; i++) {
		j = i - 1;
		current = eigval[i];
		currentVec = eigvec[i];
		past = eigval[j];
		pastVec = eigvec[j];
		while (current < past) {
			eigval[j+1] = past;
			eigvec[j+1] = pastVec;
			j--;
			if (j < 0) {
				break;
			}
			past = eigval[j];
			pastVec = eigvec[j];
		}
		eigval[j+1] = current;
		eigvec[j+1] = currentVec;
	}
}

// transpose n x n matrix
double **transpose(double **A, int n) {
	double **transA = new double*[n];
	for (int i = 0; i < n; i++) {
		transA[i] = new double[n];
		for (int j = 0; j < n; j++) {
			transA[i][j] = A[j][i];
		}
	}
	return transA;
}

double analyticConvergenceRate(int n, double eps, double sumOff) {
	double N = (1.0 - 2.0/(n*n - n));
	return ((log(eps) - log(sumOff))/log(N));
}

void normalize(double *v, double *u, int n) {
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += v[i];
	}
	for (int j = 0; j < n; j++) {
		u[j] = v[j]/sum;
	}
}

void extractEigenVec(double **A, double *u, int index, int n) {
	for (int i = 0; i < n; i++) {
		u[i] = A[index][i];
	}
}