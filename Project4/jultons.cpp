#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <random>

template <typename myType>
myType *createArray(myType value, int n);

template <typename myType>
myType **createMatrix(myType value, int row, int col);

template <typename myType>
myType arrayToFile(myType *v , int n, std::string filename);

int **createRandomSpinMatrix(int row, int col);

void mcmc(double itemp, double ftemp, double tempStep, int N, int mcs);

void metropolis(int N, int &counter, int **spinMatrix, double &E, double &M, double *w);

void output(double values, double accConf, double temp, int mcs, int N);

int periodicBoundary(int i, int limit, int add);

double ranIO();


int main() {

    return 0;
}

void mcmc(double itemp, double ftemp, double tempStep, int N, int mcs) {
    double E, M, *values, *w;
    int **spinMatrix, counter, *accConf;
    values = createArray<double>(0, 5);
    w = createArray<double>(0, 17);
    accConf = createArray<int>(0, mcs);
    spinMatrix = createMatrix<int>(1, N, N);
    for(double temp = itemp; temp < ftemp; temp += tempStep) {
        M = N*N;
        E = -2*N*N;
        for(int de =-8; de <= 8; de+=4) {
            w[de+8] = exp(-de/temp);
        }
        for(int cycles = 0; cycles < mcs; cycles++) {
            counter = 0;
            metropolis(N, counter, spinMatrix, E, M, w);
            values[0] += E;
            values[1] += E*E;
            values[2] += M;
            values[3] += M*M;
            values[4] += fabs(M);
            accConf[cycles] = counter;
        }
    }
}

void output(double *values, double *accConf, double temp, int mcs, int N) {
    double norm = 1/((double)(mcs));
    int nn = N/N;
    double Emean = values[0]*norm;
    double E2mean = values[1]*norm;
    double Mmean = values[2]*norm;
    double M2mean = values[3]*norm;
    double Mabsmean = values[4]*norm;

    double Evar = (E2mean - Emean*Emean)/nn;
    double Mvar = (M2mean - Mmean*Mmean)/nn;
    double M2var = (M2mean - Mabsmean*Mabsmean)/nn;
}

void metropolis(int N, int &counter, int **spinMatrix, double &E, double &M, double *w) {
    for(int y = 0; y < N; y++) {
        for(int x = 0; x < N; x++) {
            int ix = (int) (ranIO()*(double)N);
            int iy = (int) (ranIO()*(double)N);
            int deltaE = 2*spinMatrix[iy][ix]*
            (spinMatrix[iy][periodicBoundary(ix, N, -1)] +
            spinMatrix[periodicBoundary(iy, N, -1)][ix] +
            spinMatrix[iy][periodicBoundary(ix, N, 1)] +
            spinMatrix[periodicBoundary(iy, N, 1)][ix]);
            if(ranIO() <= w[deltaE + 8]) {
                spinMatrix[iy][ix] *= -1;
                M += (double) 2*spinMatrix[iy][ix];
                E += (double) deltaE;
                counter +=1;
            }
        }
    }
}

int periodicBoundary(int i, int limit, int add) {
    return (i+limit+add) % (limit);
}

template <typename myType>
myType **createMatrix(myType value, int row, int col) {
    myType **mat;
    mat = new myType*[row];
    for(int i = 0; i < row; i ++) {
        mat[i] = new myType[col];
        for(int j = 0; j < col; j++) {
            mat[i][j] = value;
        }
    }
    return mat;
}

int **createRandomSpinMatrix(int row, int col) {
    int **mat, value;
    mat = new int*[row];
    for(int i = 0; i < row; i ++) {
        mat[i] = new int[col];
        for(int j = 0; j < col; j++) {
            value = rand()%2;
            if (value == 0) mat[i][j] = -1;
            else mat[i][j] = value;
        }
    }
    return mat;
}

template <typename myType>
myType *createArray(myType value, int n) {
    myType *vec;
    vec = new myType[n];
    for(int i = 0; i < n; i++) {
        vec[i] = value;
    }
    return vec;
}

template <typename myType>
myType arrayToFile(myType *v , int n, std::string filename) {
	std::ofstream myfile(filename + ".txt");
	if (myfile.is_open()) {
		myfile << n << "\n";
		for (int i = 0; i < n; i++) {
			myfile << v[i] << "\n";
		}
	}
}

double ranIO() {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
    return RandomNumberGenerator(gen);
}
