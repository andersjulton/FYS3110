#include <string>
using namespace std;

double* createVector(double value, int n);

double maxError(double *expected, double *computed, int n);

double maxEpsilon(double *expected, double *computed, int n);

void printError(double *u, double *v, int n);

void arrayToFile(double *v , int n, string filename);