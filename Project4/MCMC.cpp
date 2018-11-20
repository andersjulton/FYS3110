#include <fstream>
#include <iostream>
#include <mpi.h>
#include "MCMC.h"
#include "utils.h"

// Function for handling the boundary conditions
inline int periodicBoundary(int i, int limit, int add) {
    return (i+limit+add) % (limit);
}


// Function used to evaluate implementation and test validity
void mcmc(double itemp, double ftemp, double tempStep, int L, int mcs, bool random, int preCycles, std::string filename, long idumshift) {
	double E, M, values[5], totalValues[5], w[17], **outputValues;
	int counter, **spinMatrix, j, outputInterval = mcs/10, counterNorm;
	long idum = -1 - idumshift;

	spinMatrix = (int**)createPointerMatrix(L, L, sizeof(int));
	outputValues = (double**)createPointerMatrix(7, outputInterval, sizeof(double));

	for (double temp = itemp; temp < ftemp; temp += tempStep) {
		j = 0;
		if (temp == itemp) {
			init(idum, w, values, totalValues, spinMatrix, L, mcs, temp, E, M, random, preCycles);
		}
		else {
			reachEquib(0, L, spinMatrix, idum, w, E, M);
		}
		for (int cycles = 0; cycles < mcs; cycles++) {
			counter = 0;
			metropolis(L, idum, counter, spinMatrix, E, M, w);
			// Update expectation values
			values[0] += E;
			values[1] += E*E;
			values[2] += M;
			values[3] += M*M;
			values[4] += fabs(M);
			// Extract normalized expectation values at a set interval
			if ((cycles) % 10 == 0) {
				counterNorm = cycles + 1;
				outputValues[0][j] = values[0]/counterNorm;	
				outputValues[1][j] = values[4]/counterNorm;	
				outputValues[2][j] = (values[1]/counterNorm - values[0]*values[0]/counterNorm/counterNorm);	
				outputValues[3][j] = (values[3]/counterNorm - values[4]*values[4]/counterNorm/counterNorm);	
				outputValues[4][j] = counter;
				outputValues[5][j] = E;							
				outputValues[6][j] = fabs(M);
				j++;
			}
		}
		output(values, outputValues, temp, mcs, L, outputInterval, filename);	
	}
	destroyPointerMatrix((void**)spinMatrix);
	destroyPointerMatrix((void**)outputValues);
}

// Function for simulating larger lattice sizes
void mcmcPara(double itemp, double ftemp, double tempStep, int L, int mcs, bool random, int preCycles, std::string filename) {
	double E, M, values[5], totalValues[5], w[17], *eBuffer, *eFull, **output;
	int counter, **spinMatrix, gsize, myRank, j, k = 0;
	int fint = (ftemp - itemp)/tempStep + 1.5;
	double norm = (1.0/(double(mcs)));

	MPI_Comm_size(MPI_COMM_WORLD, &gsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	if (myRank == 0) {
		std::cout << "Beginning calculations with " << gsize << " processes.\n";
		std::cout << "Lattice size: " << L << "x" << L << "\n";
		std::cout << mcs << " monte carlo cycles\n";
	}

	// Divide total number of MC cycles into intervals for each process 
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
		doubleArrayToBinary(eFull, mcs, filename+"E");
		writeMatrixDim(7, fint, filename);
		doubleMatrixToBinary(output, 7, fint, filename);
	}
    destroyPointerMatrix((void**)spinMatrix);
	destroyPointerMatrix((void**)output);
	delete[] eBuffer;
	delete[] eFull;
}

// Function for performing a Monte Carlo step with the Metropolis test
void metropolis(int L, long &idum, int &counter, int **spinMatrix, double &E, double &M, double *w) {
    for(int i = 0; i < L*L; i++) {
		// Choose random spin
		int ranj = (int) (ran2(&idum)*(double)L);
        int rani = (int) (ran2(&idum)*(double)L);

		// Calculate energy difference when flipping one spin
        int deltaE = 2*spinMatrix[rani][ranj]*
        (spinMatrix[rani][periodicBoundary(ranj, L, -1)] +
			spinMatrix[periodicBoundary(rani, L, -1)][ranj] +
			spinMatrix[rani][periodicBoundary(ranj, L, 1)] +
			spinMatrix[periodicBoundary(rani, L, 1)][ranj]);

		// Flip spin and update observables if Metropolis test is passed.
		// For dE < 0, exp(-dE/temp) > 1, so test will always be passed for these values
		if(ran2(&idum) <= w[deltaE + 8]) {
			spinMatrix[rani][ranj] *= -1;
            M += (double) 2*spinMatrix[rani][ranj];
            E += (double) deltaE;
            counter++; // For counting accepted configurations
        }
    }
}

// Writes expectation values to binary file
void output(double *totalValues, double **outputValues, double temp, int mcs, int L, int outputInterval, std::string filename) {
	// Expectation values are written to file for plotting in python.
	valueInfoToTXT(7, outputInterval, temp, filename);
	doubleMatrixToBinary(outputValues, 7, outputInterval, filename);
}

// Initializing various value arrays and the spin matrix
void init(long &idum, double *w, double *values, double *totalValues, int **spinMatrix, 
	int L, int mcs, double temp, double &E, double &M, bool random, int preCycles) {
    E = M = 0;
  
    for(int i = 0; i < 17; i++) {w[i] = 0;}
    for(int i = 0; i < 5; i++) {
		values[i] = 0;
		totalValues[i] = 0;
	}
	// Set up the array for holding the possible Boltzmann factors
    for(int de =-8; de <= 8; de+=4) {w[de+8] = exp(-de/temp);}

	// For random configuration
    if (random) {
        for(int i = 0; i < L; i++) {
            for(int j = 0; j < L; j++) {
                if(ran2(&idum)>=0.5) {spinMatrix[i][j] = 1;}
                else {spinMatrix[i][j] = -1;}
            }
        }
		reachEquib(preCycles, L, spinMatrix, idum, w, E, M);
    }

	// For ordered (all spins up) configuration
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

// Runs a set number of monte carlo cycles before measurements
void reachEquib(int preCycles, int L, int **spinMatrix, long &idum, double *w, double &E, double &M) {
	for (int i = 0; i < preCycles; i++) {
		for (int j = 0; j < L*L; j++) {
			int ranj = (int)(ran2(&idum)*(double)L);
			int rani = (int)(ran2(&idum)*(double)L);
			int deltaE = 2*spinMatrix[rani][ranj]*
				(spinMatrix[rani][periodicBoundary(ranj, L, -1)] +
					spinMatrix[periodicBoundary(rani, L, -1)][ranj] +
					spinMatrix[rani][periodicBoundary(ranj, L, 1)] +
					spinMatrix[periodicBoundary(rani, L, 1)][ranj]);
			if (ran2(&idum) <= w[deltaE + 8]) {
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


// Function for calculating energy probabilites. 
void getProb(double *a, int mcs, std::string filename) {
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
	
	valueInfoToTXT(2, k, 0, filename+"_prob");
	doubleMatrixToBinary(output, 2, k, filename+"_prob");

	delete[] seen;
	deleteMatrix(output, 2);
	deleteMatrix(buffer, 2);
}