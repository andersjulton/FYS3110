#pragma once
#include <string>

void init(long &idum, double *w, double *values, double *totaValues, int **spinMatrix,
	int L, int mcs, double temp, double &E, double &M, bool random, int preCycles);

void mcmcPara(double itemp, double ftemp, double tempStep, int L, int mcs, bool random, int preCycles, std::string filename);

void mcmc(double itemp, double ftemp, double tempStep, int L, int mcs, bool random, int preCycles, std::string filename, long idumshift);

void metropolis(int L, long &idum, int &counter, int **spinMatrix, double &E, double &M, double *w);

void output(double *values, double **outputValues, double temp, int mcs, int L, int outputInterval, std::string filename);

//void outputPara(double **totalValues, double *eFull,  double temp, int mcs, int L, int k, std::string filename);

void valueInfoToTXT(int n, int m, double temp, std::string filename);

void reachEquib(int cycles, int L, int **spinMatrix, long &idum, double *w, double &E, double &M);

void getProb(double *a, int n, std::string filename);

template <typename myType>
myType *createArray(myType value, int n);