#include <string>

using namespace std;


double* diagToVector(double **A, int n);

double* createVector(double value, int n);

double **createMatrix(int m, int n);

double **createDiaMatrix(double d, int n);

double **createTriDiaMatrix(double off_value, double d_value, int n);

void deleteMatrix(double **mat, int n);

double maxError(double *expected, double *computed, int n);

double maxEpsilon(double *expected, double *computed, int n);

void printError(double *u, double *v, int n);

void arrayToFile(double *v , int n, string filename);
