#include <iostream>
#include <string>
#include <mpi.h>
#include "analyticValues.h"
#include "MCMC.h"
	
void b();
void c();
void d();
void e();

int main(int argc, char* argv[]) {

	//b();
	//c();

	//MPI_Init(&argc, &argv);
	//d();
	//e();
	//MPI_Finalize();

	return 0;
}

void b() {
	double itemp = 1.0;
	double ftemp = 1.1;
	double tempStep = 0.1;
	int L = 2;
	int mcs = int(1e7);
	int preCycles = 0;
	bool random = true;
	long idumshift = 5;
	for (int i = 0; i < 1; i++) {
		idumshift = i;
		std::string filename = "ex_b"+std::to_string(i);
		mcmc(itemp, ftemp, tempStep, L, mcs, random, preCycles, filename, idumshift);
	}
}

void c() {
	double itemp = 1.0;
	double ftemp = 1.1;
	double tempStep = 0.1;
	int L = 20;
	int mcs = int(1e5);
	int preCycles = 0;
	bool random = true;
	long idumshift = 5;

	std::string filename = "ex_c_1_unord";
	mcmc(itemp, ftemp, tempStep, L, mcs, random, preCycles, filename, idumshift);

	random = false;
	filename = "ex_c_1_ord";
	mcmc(itemp, ftemp, tempStep, L, mcs, random, preCycles, filename, idumshift);

	itemp = 2.4;
	ftemp = 2.5;
	random = true;
	filename = "ex_c_2.4_unord";
	mcmc(itemp, ftemp, tempStep, L, mcs, random, preCycles, filename, idumshift);

	random = false;
	filename = "ex_c_2.4_ord";
	mcmc(itemp, ftemp, tempStep, L, mcs, random, preCycles, filename, idumshift);
}

void d() {

	double itemp = 1.0;
	double ftemp = 1.1;
	double tempStep = 0.1;
	int L = 40;
	int mcs = (int) 1e6;
	int preCycles = 5000;
	bool random = true;

	std::string filename = "ex_d_1_REV";
	mcmcPara(itemp, ftemp, tempStep, L, mcs, random, preCycles, filename);

	itemp = 2.4;
	ftemp = 2.5;
	filename = "ex_d_2.4_REV";
	mcmcPara(itemp, ftemp, tempStep, L, mcs, random, preCycles, filename);
}

void e() {

	double itemp = 2.25;
	double ftemp = 2.29;
	double tempStep = 0.0015;
	int L = 40;
	int mcs = (int) 1e6;
	int preCycles = 5000;
	bool random = true;

	std::string filename = "ex_e_40x40_2.25-2.29";
	mcmcPara(itemp, ftemp, tempStep, L, mcs, random, preCycles, filename);

	L = 60;
	filename = "ex_e_60x60_2.25-2.29";
	mcmcPara(itemp, ftemp, tempStep, L, mcs, random, preCycles, filename);

	L = 80;
	filename = "ex_e_80x80_2.25-2.29";
	mcmcPara(itemp, ftemp, tempStep, L, mcs, random, preCycles, filename);

	L = 100;
	filename = "ex_e_100x100_2.25-2.29";
	mcmcPara(itemp, ftemp, tempStep, L, mcs, random, preCycles, filename);

	L = 140;
	filename = "ex_e_140x140_2.5-2.29";
	mcmcPara(itemp, ftemp, tempStep, L, mcs, random, preCycles, filename);
}