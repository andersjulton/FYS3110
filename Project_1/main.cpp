#include <fstream>
#include "tasks.h"

int main(int argc, char* argv[]) {
	// Reading in dimensions and values for vectors
	//int n = atoi(argv[1]);
	double a_value = atoi(argv[1]);
	double b_value = atoi(argv[2]);
	double c_value = atoi(argv[3]);

	partB(a_value, b_value, c_value);

	//partC(a_value, b_value, c_value);

	//partD();

	//partE();

	return 0;
}