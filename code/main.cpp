#include <fstream>
#include "func.h"
#include <iostream>
#include <string>

void partB(int a_value, int b_value, int c_value);
void partC(int a_value, int b_value, int c_value);
void partD();
void partE();

int main(int argc, char* argv[]) {
	// Reading in dimensions and values for vectors;
	double a_value = atoi(argv[1]);
	double b_value = atoi(argv[2]);
	double c_value = atoi(argv[3]);

	printf("\n");
	std::string choice;
	std::string info = "Write any of the following commands:\n \
	B    : run only part B of the assignment\n \
	C    : run only part C of the assignment\n \
	D    : run only part D of the assignment\n \
	E    : run only part E of the assignment\n \
	all  : run all parts of the assignment\n \
	a    : change the value of variable a \n \
	b    : change the value of variable b \n \
	c    : change the value of variable c \n \
	info : write out this message\n \
	q    : quit program\n\n"; 

	std::cout << info;
	std::string number;

	int pigsThatFly = 0;

	while (pigsThatFly == 0) {
		std::cout << ">> ";
		getline(std::cin, choice);
		if (choice.compare("B") == 0) {
			partB(a_value, b_value, c_value);
		} else if (choice.compare("C") == 0) {
			partC(a_value, b_value, c_value);
		} else if (choice.compare("D") == 0) {
			partD();
		} else if (choice.compare("E") == 0) {
			partE();
		} else if (choice.compare("all") == 0) {
			partB(a_value, b_value, c_value);
			partC(a_value, b_value, c_value);
			partD();
			partE();
		} else if (choice.compare("q") == 0) {
			pigsThatFly = 1;
			printf("Bye!\n");
		} else if ((choice.compare("info") == 0) || (choice.compare("help") == 0)) {
			std::cout << info;
		} else if (choice.compare("a") == 0) {
			std::cout << "a = ";
			getline(std::cin, number);
			a_value = std::stod(number);
			printf("%.5f\n", a_value);
			printf("\n");
		} else if (choice.compare("b") == 0){
			std::cout << "b = ";
			getline(std::cin, number);
			b_value = std::stod(number);
			printf("%.5f\n", b_value);
			printf("\n");
		} else if (choice.compare("c") == 0) {
			std::cout << "c = ";
			getline(std::cin, number);
			c_value = std::stod(number);
			printf("%.5f\n", c_value);
			printf("\n");
		} else {	
			std::cout << "not a valid choice\n\n";
		}
	}
	return 0;
}

void partB(int a_value, int b_value, int c_value) {
	printf("Project 1, part B\n");
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
		printf("File n_%d.txt has been written.\n", n);

		delete[] a;
		delete[] b;
		delete[] c;
		delete[] f;
		delete[] v;

		n = n*10;
	}
	printf("\n\n");
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