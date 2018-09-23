#include <fstream>
#include "func.h"
#include "utils.h"
#include "testFunc.h"
#include "armadillo"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <string>

using namespace std;
using namespace arma;

int main() {

	return 0;
} 

void partC(int a_value, int b_value, int c_value) {
	printf("Project 1, part C\n");
	int n = 10;
	for (int i = 0; i < 8; i++) {
		compareTime(a_value, b_value, c_value, n);
		n = n*10;
	}
	printf("\n");
}

void partE() {
	printf("Project 1, part E\n");
	int n = 10;
	for (int i = 0; i < 4; i++) {
		compareTimeArmadillo(n);
		n = n*10;
	}
	//compareTimeArmadillo(n);
	printf("\n");
}

void partD() {
	printf("Project 1, part D\n\n");
	double error;
	double h;
	int n = 10;
	printf("         n |       log10(h)       |       max error\n");
	printf("   ---------------------------------------------------\n");
	for (int i = 0; i < 7; i++) {
		h = 1.0/(n + 1.0);
		error = maxErrorDiaSolver(n);
		printf("%10d | %20.15f | %18.15f\n", n, log10(h), error);
		n = n*10;
	}
	printf("\n");
}

