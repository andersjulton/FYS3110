#include "pch.h"  
#include <cmath>
#include <stdio.h>
#include <fstream>
#include "utils.h"
#include "thomas.h"
#include "tests.h"
#include "tasks.h"

int main(int argc, char* argv[]) {
	// Reading in dimensons and values for vectors
	//int n = atoi(argv[1]);
	double a_value = atoi(argv[2]);
	double b_value = atoi(argv[3]);
	double c_value = atoi(argv[4]);

	partB(a_value, b_value, c_value);

	partC(a_value, b_value, c_value);

	partD();

	partE();

	return 0;
}