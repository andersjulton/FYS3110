#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <random>


template <typename myType>
myType *createArray(myType value, int n);

template <typename myType>
myType **createMatrix(int row, int col, myType value);

template <typename myType>
myType arrayToFile(myType *v , int n, std::string filename);

void **createPointerMatrix(int row, int col, int num_bytes);

void init(int *accConf, double *w, double *values, int **spinMatrix, int N, int mcs, double temp, double &E, double &M, bool random);

void mcmc(double itemp, double ftemp, double tempStep, int N, int mcs, bool random);

void metropolis(int N, int &counter, int **spinMatrix, double &E, double &M, double *w);

void output(double *values, int *accConf, double temp, int mcs, int N);

int periodicBoundary(int i, int limit, int add);



int main() {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

    double itemp = 1;
    double ftemp = 2;
    double tempStep = 1;
    int N = 4;
    int mcs = 8;
    bool random = true;

    mcmc(itemp, ftemp, tempStep, N, mcs, random);
    return 0;
}

void mcmc(double itemp, double ftemp, double tempStep, int N, int mcs, bool random) {
    double E, M, values[5], w[17];
    int counter, accConf[mcs], **spinMatrix;
    spinMatrix = (int**) createPointerMatrix(N, N, sizeof(int));
    for(double temp = itemp; temp < ftemp; temp += tempStep) {
        init(accConf, w, values, spinMatrix, N, mcs, temp, E, M, random);
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
        output(values, accConf, temp, mcs, N);
    }
}

void metropolis(int N, int &counter, int **spinMatrix, double &E, double &M, double *w) {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
    for(int y = 0; y < N; y++) {
        for(int x = 0; x < N; x++) {
            int ix = (int) (RandomNumberGenerator(gen)*(double)N);
            int iy = (int) (RandomNumberGenerator(gen)*(double)N);
            int deltaE = 2*spinMatrix[iy][ix]*
            (spinMatrix[iy][periodicBoundary(ix, N, -1)] +
            spinMatrix[periodicBoundary(iy, N, -1)][ix] +
            spinMatrix[iy][periodicBoundary(ix, N, 1)] +
            spinMatrix[periodicBoundary(iy, N, 1)][ix]);
            if(RandomNumberGenerator(gen) <= w[deltaE + 8]) {
                spinMatrix[iy][ix] *= -1;
                M += (double) 2*spinMatrix[iy][ix];
                E += (double) deltaE;
                counter += 1;
            }
        }
    }
}

void output(double *values, int *accConf, double temp, int mcs, int N) {
    double norm = 1/((double)(mcs));
    double Emean = values[0]*norm;
    double E2mean = values[1]*norm;
    double Mmean = values[2]*norm;
    double M2mean = values[3]*norm;
    double Mabsmean = values[4]*norm;

    double Evar = (E2mean - Emean*Emean)/N/N;
    double Mvar = (M2mean - Mmean*Mmean)/N/N;
    double M2var = (M2mean - Mabsmean*Mabsmean)/N/N;
    std::cout << Evar << '\n';

}

void init(int *accConf, double *w, double *values, int **spinMatrix, int N, int mcs, double temp, double &E, double &M, bool random) {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
    E = M = 0;
    for(int i = 0; i < mcs; i++) {accConf[i] = 0;}
    for(int i = 0; i < 17; i++) {w[i] = 0;}
    for(int i = 0; i < 5; i++) {values[i] = 0;}
    for(int de =-8; de <= 8; de+=4) {w[de+8] = exp(-de/temp);}
    if (random) {
        int value;
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                value = rand()%2;
                if (value == 0) {spinMatrix[i][j] = -1;}
                else {spinMatrix[i][j] = value;}
                M += (double) spinMatrix[i][j];
            }
        }
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                E -= (double) spinMatrix[i][j]*
                (spinMatrix[periodicBoundary(i, N, -1)][j] +
                spinMatrix[i][periodicBoundary(j, N, -1)]);
            }
        }
    }
    else {
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                spinMatrix[i][j] = 1;
            }
        }
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                E -= (double) spinMatrix[i][j]*
                (spinMatrix[periodicBoundary(j, N, -1)][i] +
                spinMatrix[j][periodicBoundary(i, N, -1)]);
            }
        }
        M = N*N;
        E = -2*N*N;
    }
}

int periodicBoundary(int i, int limit, int add) {
    return (i+limit+add) % (limit);
}



int **createRandomSpinMatrix(int row, int col, int &M) {
    int **mat, value;
    mat = new int*[row];
    for(int i = 0; i < row; i ++) {
        mat[i] = new int[col];
        for(int j = 0; j < col; j++) {
            value = rand()%2;
            if (value == 0) mat[i][j] = -1;
            else mat[i][j] = value;
            M += mat[i][j];
        }
    }
    return mat;
}

template <typename myType>
myType **createMatrix(int row, int col, myType value) {
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

void **createPointerMatrix(int row, int col, int num_bytes) {
    int i, num;
    char **pointer, *ptr;

    pointer = new(std::nothrow) char*[row];
    if(!pointer) {
        std::cout << "Exception handling: Memory allocation failed";
        std::cout << " for "<< row << "row addresses !" << '\n';
        return NULL;
    }
    i = (row*col*num_bytes)/sizeof(char);
    pointer[0] = new(std::nothrow) char[i];
    if(!pointer[0]) {
        std::cout << "Exception handling: Memory allocation failed";
        std::cout << " for address to " << i << " characters !" << '\n';
        return NULL;
    }
    ptr = pointer[0];
    num = col*num_bytes;
    for(i = 0; i < row; i++, ptr += num )   {
        pointer[i] = ptr;
    }
    return (void **)pointer;
}
