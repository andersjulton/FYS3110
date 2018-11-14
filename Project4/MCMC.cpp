#include <fstream>
#include <iostream>
#include <mpi.h>
#include "MCMC.h"
#include "utils.h"


inline int periodicBoundary(int i, int limit, int add) {
    return (i+limit+add) % (limit);
}

void mcmc(double itemp, double ftemp, double tempStep, int L, int mcs, bool random, int preCycles, std::string filename, long idumshift) {
	double E, M, values[5], totalValues[5], w[17], **outputValues;
	int counter, **spinMatrix, j, outputInterval = mcs/10, counterNorm;
	long idum = -1 - idumshift;

	spinMatrix = (int**)createPointerMatrix(L, L, sizeof(int));
	outputValues = (double**)createPointerMatrix(6, outputInterval, sizeof(double));

	for (double temp = itemp; temp < ftemp; temp += tempStep) {
		j = 0;
		init(idum, w, values, totalValues, spinMatrix, L, mcs, temp, E, M, random, preCycles);
				
		for (int cycles = 0; cycles < mcs; cycles++) {
			counter = 0;
			metropolis(L, idum, counter, spinMatrix, E, M, w);
			values[0] += E;
			values[1] += E*E;
			values[2] += M;
			values[3] += M*M;
			values[4] += fabs(M);
			if ((cycles) % 10 == 0) {
				counterNorm = cycles + 1;
				outputValues[0][j] = values[0]/counterNorm;	// mean energy
				outputValues[1][j] = values[4]/counterNorm;	// mean abs magn
				outputValues[2][j] = (values[1]/counterNorm - values[0]*values[0]/counterNorm/counterNorm);	// CV
				outputValues[3][j] = (values[3]/counterNorm - values[4]*values[4]/counterNorm/counterNorm);	// Susc
				outputValues[4][j] = counter;					// accepted configs per MC cycle
				outputValues[5][j] = E;							// for prob measurements
				j++;
			}
		}
		output(values, outputValues, temp, mcs, L, outputInterval, filename);	
	}
	destroyPointerMatrix((void**)spinMatrix);
	destroyPointerMatrix((void**)outputValues);
}

void mcmcPara(double itemp, double ftemp, double tempStep, int L, int mcs, bool random, int preCycles, std::string filename) {
	double E, M, values[5], totalValues[5], w[17], *eBuffer, *eFull, **output;
	int counter, **spinMatrix, gsize, myRank, j, k = 0;
	int fint = (ftemp - itemp)/tempStep + 1.5;
	double norm = (1.0/(double(mcs)));

	MPI_Comm_size(MPI_COMM_WORLD, &gsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	if (myRank == 0) {
		std::cout << "Beginning calculations with " << gsize << " processors.\n";
		std::cout << "Lattice size: " << L << "x" << L << "\n";
		std::cout << mcs << " monte carlo cycles\n";
	}

	// Divide total number of MC cycles into intervals for each processor 
	int intervals = mcs/gsize;
	int istart = myRank*intervals;
	int iend = (myRank+1)*intervals;
	if ((myRank == gsize - 1) && (iend < mcs)) iend = mcs;

	// Arrays for storing energies for probability measurements
	eBuffer = new double[intervals];
	eFull = new double[mcs];

	MPI_Bcast(&L, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&itemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ftemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&tempStep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    spinMatrix = (int**) createPointerMatrix(L, L, sizeof(int));
	output = (double**)createPointerMatrix(7, fint, sizeof(double));

	long idum = -3 - myRank;
	double iTime, fTime, totTime;
	iTime = MPI_Wtime();

    for(double temp = itemp; temp < ftemp; temp += tempStep) {
		j = 0;
        init(idum, w, values, totalValues, spinMatrix, L, mcs, temp, E, M, random, preCycles);
		for (int cycles = istart; cycles < iend; cycles++) {
			counter = 0;
			metropolis(L, idum, counter, spinMatrix, E, M, w);
			values[0] += E;
			values[1] += E*E;
			values[2] += M;
			values[3] += M*M;
			values[4] += fabs(M);
			eBuffer[j] = E;
			j++;
		}

		for (int i = 0; i < 5; i++) {
			MPI_Reduce(&values[i], &totalValues[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		}

		// Mean values per spin
		if (myRank == 0) {
			for (int i = 0; i < 5; i++) {
				output[i][k] = totalValues[i]*norm/L/L;
			}
			output[5][k] = (totalValues[1]*norm - totalValues[0]*totalValues[0]*norm*norm)/L/L;
			output[6][k] = (totalValues[3]*norm - totalValues[4]*totalValues[4]*norm*norm)/L/L;
		}
		
		MPI_Gather(eBuffer, intervals, MPI_DOUBLE, eFull, intervals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		k++;
    }
	fTime = MPI_Wtime();
	totTime = fTime - iTime;
	if (myRank == 0) {
		
		std::cout << "Elapsed time = " << totTime << "\n";
		writeMatrixDim(5, fint, filename);
		doubleMatrixToBinary(output, 7, fint, filename);
	}
    destroyPointerMatrix((void**)spinMatrix);
	destroyPointerMatrix((void**)output);
	delete[] eBuffer;
	delete[] eFull;
}

void metropolis(int L, long &idum, int &counter, int **spinMatrix, double &E, double &M, double *w) {
    for(int i = 0; i < L*L; i++) {
		int ranj = (int) (ran1(&idum)*(double)L);
        int rani = (int) (ran1(&idum)*(double)L);
        int deltaE = 2*spinMatrix[rani][ranj]*
        (spinMatrix[rani][periodicBoundary(ranj, L, -1)] +
			spinMatrix[periodicBoundary(rani, L, -1)][ranj] +
			spinMatrix[rani][periodicBoundary(ranj, L, 1)] +
			spinMatrix[periodicBoundary(rani, L, 1)][ranj]);
		if(ran1(&idum) <= w[deltaE + 8]) {
			spinMatrix[rani][ranj] *= -1;
            M += (double) 2*spinMatrix[rani][ranj];
            E += (double) deltaE;
            counter++;
        }
    }
}

/*


void outputPara(double **totalValues, double *eFull, double temp, int mcs, int L, int k, std::string filename) {
	double norm = (1.0/(double(mcs)));
	double Emean = totalValues[0]*norm;
	double E2mean = totalValues[1]*norm;
	double Mmean = totalValues[2]*norm;
	double M2mean = totalValues[3]*norm;
	double Mabsmean = totalValues[4]*norm;
	double Evar = (E2mean - Emean*Emean);
	double M2var = (M2mean - Mabsmean*Mabsmean);
	double CV = Evar/temp/temp;
	double susc = M2var/temp;
	
	//getProb(eFull, mcs, temp, filename);
	outputValues[0][k] = Emean;
	outputValues[1][k] = Mabsmean;
	outputValues[2][k] = CV;
	outputValues[3][k] = susc;
}
*/
void output(double *totalValues, double **outputValues, double temp, int mcs, int L, int outputInterval, std::string filename) {
    double norm = (1.0/(double(mcs)));
    double Emean = totalValues[0]*norm;
    double E2mean = totalValues[1]*norm;
    double Mmean = totalValues[2]*norm;
    double M2mean = totalValues[3]*norm;
    double Mabsmean = totalValues[4]*norm;

    double Evar = (E2mean - Emean*Emean)/temp/temp;
    double M2var = (M2mean - Mabsmean*Mabsmean)/temp;
	/*
	std::cout << "\nNumerical values for temperature = " << temp << " and " << mcs << " monte carlo cycles\n";
	std::cout << "Mean energy = " << Emean << '\n';
	std::cout << "Mean energy squared = " << E2mean << '\n';
	std::cout << "Mean absolute magnetization = " << Mabsmean << '\n';
	std::cout << "Mean magnetization squared = " << M2mean << '\n';
	std::cout << "Specific heat = " << Evar << '\n';
	std::cout << "Magnetic susceptibility = " << M2var << '\n';*/
	std::string str = std::to_string(temp);
	str.erase(str.find_last_not_of('0') + 1, std::string::npos);
	valueInfoToTXT(6, outputInterval, temp, filename);
	doubleMatrixToBinary(outputValues, 6, outputInterval, filename);
}

void init(long &idum, double *w, double *values, double *totalValues, int **spinMatrix, 
	int L, int mcs, double temp, double &E, double &M, bool random, int preCycles) {
    E = M = 0;
  
    for(int i = 0; i < 17; i++) {w[i] = 0;}
    for(int i = 0; i < 5; i++) {
		values[i] = 0;
		totalValues[i] = 0;
	}
    for(int de =-8; de <= 8; de+=4) {w[de+8] = exp(-de/temp);}

    if (random) {
        for(int i = 0; i < L; i++) {
            for(int j = 0; j < L; j++) {
                if(ran1(&idum)>=0.5) {spinMatrix[i][j] = 1;}
                else {spinMatrix[i][j] = -1;}
            }
        }
		reachEquib(preCycles, L, spinMatrix, idum, w, E, M);
    }
    else {
        for(int i = 0; i < L; i++) {
            for(int j = 0; j < L; j++) {
                spinMatrix[i][j] = 1;
            }
        }
		E = -2*L*L;
		M = L*L;
    }
}

void reachEquib(int preCycles, int L, int **spinMatrix, long &idum, double *w, double &E, double &M) {
	// Runs a set number of monte carlo cycles before measurements
	for (int i = 0; i < preCycles; i++) {
		for (int j = 0; j < L*L; j++) {
			int ranj = (int)(ran1(&idum)*(double)L);
			int rani = (int)(ran1(&idum)*(double)L);
			int deltaE = 2*spinMatrix[rani][ranj]*
				(spinMatrix[rani][periodicBoundary(ranj, L, -1)] +
					spinMatrix[periodicBoundary(rani, L, -1)][ranj] +
					spinMatrix[rani][periodicBoundary(ranj, L, 1)] +
					spinMatrix[periodicBoundary(rani, L, 1)][ranj]);
			if (ran1(&idum) <= w[deltaE + 8]) {
				spinMatrix[rani][ranj] *= -1;
			}
		}
	}
	// Retrieve starting energy and magnetisation
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			M += (double)spinMatrix[i][j];
			E -= (double)spinMatrix[i][j]*
				(spinMatrix[i][periodicBoundary(j, L, -1)] +
					spinMatrix[periodicBoundary(i, L, -1)][j]);
		}
	}
}

template <typename myType>
myType *createArray(myType value, int n) {
	myType *vec;
	vec = new myType[n];
	for (int i = 0; i < n; i++) {
		vec[i] = value;
	}
	return vec;
}

void valueInfoToTXT(int n, int m,  double temp, std::string filename) {
	std::ofstream file(filename + "_INFO.txt");
	file << n << " " << m << " " << temp;
	file.close();
}

void getProb(double *a, int mcs, double temp, std::string filename) {
	int *seen, k = 0;
	double **buffer, **output;

	buffer = createMatrix(2, mcs);
	seen = createArray(0, mcs);

	// Counts multiplicity of distinct energies and calculates probability
	for (int i = 0; i < mcs; i++) {
		if (seen[i] == 0) {
			int count = 0;
			for (int j = i; j < mcs; j++)
				if (a[j] == a[i]) {
					count += 1;
					seen[j] = 1;
				}
			buffer[0][k] = a[i];
			buffer[1][k] = (double)count/((double)mcs);
			k++;
		}
	}
	output = createMatrix(2, k+1);
	for (int i = 0; i < k; i++) {
		output[0][i] = buffer[0][i];
		output[1][i] = buffer[1][i];
	}

	std::string str = std::to_string(temp);
	str.erase(str.find_last_not_of('0') + 1, std::string::npos);
	valueInfoToTXT(2, k, temp, filename);
	doubleMatrixToBinary(output, 2, k, filename);

	delete[] seen;
	deleteMatrix(output, 2);
	deleteMatrix(buffer, 2);
}