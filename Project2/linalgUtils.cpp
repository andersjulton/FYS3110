#include <fstream>
#include "utils.h"
#include "armadillo"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <string>

using namespace std;
using namespace arma;

/* version of insertion sort for sorting a list of eigenvalues 
their corresponding eigenvectors will be placed accordingly
eigenvec must be a list of columns where each colomn is an eigenvector*/
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
double** transpose(double **A, int n) {
	double **transA = new double*[n];
	for (int i = 0; i < n; i++) {
		transA[i] = new double[n];
		for (int j = 0; j < n; j++) {
			transA[i][j] = A[j][i];
		}
	}
	return transA;
}

// copy symmetric matrix A to arma mat
mat copySymMatrixToArma(double **A, int n) {
	mat armaA(n, n, fill::zeros);

	for (int i = 0; i < n; i++) {
		armaA(i, i) = A[i][i];
		for (int j = i+1; j < n; j++) {
			armaA(i, j) = A[i][j];
			armaA(j, i) = A[i][j];
		}
	}
	return armaA;
}









