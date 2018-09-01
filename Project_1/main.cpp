#include <fstream>
#include "tasks.h"
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
	// Reading in dimensions and values for vectors
	//int n = atoi(argv[1]);
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
		} else if (choice.compare("info") == 0) {
			std::cout << info;
		} else if (choice.compare("a") == 0) {
			std::cout << "a = ";
			getline(std::cin, number);
			a_value = std::stod(number);
			printf("%.5f\n", a_value);
			printf("\n");
		} else if (choice.compare("b") == 0) {
			std::cout << "b = ";
			getline(std::cin, number);
			b_value = std::stod(number);
			printf("%.5f\n", b_value);
			printf("\n");
		}else if (choice.compare("c") == 0) {
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