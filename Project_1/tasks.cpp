#include "tasks.h"
#include <stdio.h>
#include "utils.h"
#include "tests.h"
#include "thomas.h"
#include <math.h>

void partB(int a_value, int b_value, int c_value) {
	printf("Project 1, part B\n\n");
	double *a, *b, *c, *f, *v;
	int n = 10;

	for (int i = 0; i < 3; i++) {
		a = createVector(a_value, n);
		b = createVector(b_value, n);
		c = createVector(c_value, n);
		v = createVector(0.0, n);
		f = solutionVector(n);

		generalTriDiaSolver(a, b, c, v, f, n);
		writeToFile(v, n);

		delete[] a;
		delete[] b;
		delete[] c;
		delete[] f;
		delete[] v;

		n = n*10;
	}
}

void partC(int a_value, int b_value, int c_value) {
	printf("Project 1, part C\n");
	int n = 10;
	for(int i = 0; i < 6; i++) {
		compareTime(a_value, b_value, c_value, n);
		n = n*10;
	}
	printf("\n");
}

void partD() {
	printf("Project 1, part C\n");
	double error;
	double h;
	int n = 10;
	printf("     log10(h)                     max error\n");
	for (int i = 0; i < 7; i++) {
		h = 1.0/(n + 1.0);
		error = maxErrorDiaSolver(n);
		printf("%20.15f %25.15f\n", log10(h), error);
		n = n*10;
	}
	printf("\n");
}

void partE() {
	printf("Project 1, part E\n");
	int n = 10;
	for(int i = 0; i < 6; i++) {
		compareTimeArmadillo(n);
		n = n*10;
	}
	printf("\n");
}