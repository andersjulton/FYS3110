#include <fstream>
#include "func.h"
#include "utils.h"
#include <iostream>
#include <string>

void partB(int a_value, int b_value, int c_value);

int main(int argc, char* argv[]) {
	// Reading in dimensions and values for vectors;
	double a_value = atoi(argv[1]);
	double b_value = atoi(argv[2]);
	double c_value = atoi(argv[3]);

	printf("\n");
	std::string choice;
	std::string info = "Write any of the following commands:\n \
	B    : run only part B of the assignment\n \
	D    : run only part D of the assignment\n \
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
		} else if (choice.compare("all") == 0) {
			partB(a_value, b_value, c_value);
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
	string filename;
	for (int i = 0; i < 3; i++) {
		a = createVector(a_value, n-1);
		b = createVector(b_value, n);
		c = createVector(c_value, n-1);
		v = createVector(0.0, n);
		f = solutionVector(n);
		filename = "n_" + to_string(n);
		generalTriDiaSolver(a, b, c, v, f, n);
		arrayToFile(v, n, filename);
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
