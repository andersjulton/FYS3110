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

	b();
	//c();
	//d();
	//e();

	return 0;
}

void b() {
	double itemp = 1.0;
	double ftemp = 1.1;
	double tempStep = 0.1;
	int L = 2;
	int mcs = int(1e5);
	int preCycles = 0;
	bool random = true;
	long idumshift;
	for (int i = 0; i < 5; i++) {
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
	int mcs = 5000;
	int preCycles = 0;
	bool random = true;
	long idumshift = 0;

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
	MPI_Init(NULL, NULL);

	double itemp = 1.0;
	double ftemp = 1.1;
	double tempStep = 0.1;
	int L = 20;
	int mcs = (int) 1e5;
	int preCycles = 1000;
	bool random = true;

	std::string filename = "ex_d_1";
	mcmcPara(itemp, ftemp, tempStep, L, mcs, random, preCycles, filename);

	itemp = 2.4;
	ftemp = 2.5;
	filename = "ex_d_2.4";
	//mcmcPara(itemp, ftemp, tempStep, L, mcs, random, preCycles, filename);

	MPI_Finalize();

}

void e() {
	MPI_Init(NULL, NULL);

	double itemp = 2.2;
	double ftemp = 2.3;
	double tempStep = 0.005;
	int L = 40;
	int mcs = (int) 1e5;
	int preCycles = 1000;
	bool random = true;

	std::string filename = "ex_e_40x40_0.005";
	mcmcPara(itemp, ftemp, tempStep, L, mcs, random, preCycles, filename);

	L = 60;
	filename = "ex_e_60x60_0.005";
	//mcmcPara(itemp, ftemp, tempStep, L, mcs, random, preCycles, filename);

	L = 80;
	filename = "ex_e_80x80_0.005";
	//mcmcPara(itemp, ftemp, tempStep, L, mcs, random, preCycles, filename);

	L = 100;
	filename = "ex_e_100x100_0.005";
	//mcmcPara(itemp, ftemp, tempStep, L, mcs, random, preCycles, filename);
	MPI_Finalize();
}