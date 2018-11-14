#include <string>
#include <cmath>
#include <fstream>
#include <iostream>


//  Allocating space for a vector and fill elements with a constant value
double *createVector(double value, int n) {
	double *vector;
	vector = new double[n];
	for (int i = 0; i < n; i++) {
		vector[i] = value;
	}
	return vector;
}

// You know, like in numpy
double *linspace(double min, double max, int n) {
    double *v;
    v = new double[n];
    double step = (max - min)/(n-1);
    v[0] = min;
	for (int i = 1; i < n; i++) {
        v[i] = min + i*step;
	}
    return v;
}

int *intLinspace(int min, int max, int n) {
	int *v;
	v = new int[n];
	int step = (max - min)/(n - 1);
	v[0] = min;
	for (int i = 1; i < n; i++) {
		v[i] = min + i*step;
	}
	return v;
}

//  Allocating space for a m x n matrix and fill elements with 0
double **createMatrix(int m, int n) {
	double **mat;
	mat = new double*[m];
	for(int i = 0; i < m; i++) {
		mat[i] = new double[n];
		for(int j = 0; j < n; j++) {
			mat[i][j] = 0.0;
		}
	}
	return mat;
}

// Allocating space for a n x n diagonal matrix and fill the diagonal with d
double **createDiaMatrix(double d, int n) {
	double **mat;
	mat = createMatrix(n, n);
	// two seperate loops to avoid if-tests
	for(int i = 0; i < n; i++) {
		mat[i][i] = d;
	}
	return mat;
}

// Create a tridiagonal matrix
double **createTriDiaMatrix(double off_value, double d_value, int n) {
	double **mat;
	mat = createMatrix(n, n);
	mat[0][0] = d_value;
	mat[0][1] = off_value;
	// two seperate loops to avoid if-tests
	for (int i = 1; i < (n - 1); i++) {
		mat[i][i-1] = off_value;
		mat[i][i] = d_value;
		mat[i][i+1] = off_value;
	}
	mat[n-1][n-2] = off_value;
	mat[n-1][n-1] = d_value;
	return mat;
}

// delete 2D-array
void deleteMatrix(double **mat, int n) {
	for (int i = 0; i < n; i++) {
		delete[] mat[i];
	}
	delete[] mat;
}

// test if two arrays are equal
int testVector(double *u, double *v, int n) {
	double error;
	for (int i = 0; i < n; i++) {
		error = u[i] - v[i];
		if (fabs(error) > 1e-10) {
			printf("(%d), error = %.30f\n", i, error);
			return 1;
		}
	}
	return 0;
}

// returns the largest relative error comparing two arrays
double maxError(double *expected, double *computed, int n) {
	double max, max_i;
	max = fabs((expected[0] - computed[0])/expected[0]);
	for (int i = 1; i < n; i++) {
		max_i = fabs((expected[i] - computed[i])/expected[i]);
		if (max < max_i) {
			max = max_i;
		}
	}
	return max;
}

// returns log10 of the max relative error
double maxEpsilon(double *expected, double *computed, int n) {
	double max = maxError(expected, computed, n);
	return log10(max);
}

// write out difference between u-array and v-array
void printError(double *u, double *v, int n) {
	double error;
	printf("\n");
	for (int i = 0; i < n; i++) {
		error = u[i] - v[i];
		printf("v = %.20f, u = %.20f, error = %.30f\n", v[i], u[i], error);
	}
	printf("\n");
}

// write single array (double) to a txt-file
void doubleArrayToFile(double *v , int n, std::string filename, bool zeroPadding = false) {
	std::ofstream myfile(filename + ".txt");
	if (myfile.is_open()) {
		if (zeroPadding) {
			myfile << n+2 << "\n";
			myfile << 0.0 << "\n";
		} else {
			myfile << n << "\n";
		}
		for (int i = 0; i < n; i++) {
			myfile << v[i] << "\n";
		}
		if (zeroPadding) {
			myfile << 0.0 << "\n";
		}
	}
}

// write single array (int) to a txt-file
void intArrayToFile(int *v , int n, std::string filename, bool zeroPadding = false) {
	std::ofstream myfile(filename + ".txt");
	if (myfile.is_open()) {
		if (zeroPadding) {
			myfile << n+2 << "\n";
			myfile << 0.0 << "\n";
		} else {
			myfile << n << "\n";
		}
		for (int i = 0; i < n; i++) {
			myfile << v[i] << "\n";
		}
		if (zeroPadding) {
			myfile << 0.0 << "\n";
		}
	}
}

// write n x m matrix (double) to a txt file
void doubleMatrixToFile(double **v , int n, int m, std::string filename) {
	std::ofstream myfile(filename + ".txt");
	if (myfile.is_open()) {
		myfile << n << " " <<  m << '\n';
		for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                myfile << v[i][j] << " ";
            }
            myfile << '\n';
		}
	}
}

void doubleArrayToBinary(double *a, int n,  std::string filename) {
    std::ofstream file(filename, std::ios::binary|std::ios::out);
    file.write(reinterpret_cast<char *>(a), n*sizeof(double));
    file.close();
}

void intArrayToBinary(int *a, int n, std::string filename) {
	std::ofstream file(filename, std::ios::binary|std::ios::out);
	file.write(reinterpret_cast<char *>(a), n*sizeof(int));
	file.close();
}

// write n x m matrix (double) to a binary file
void doubleMatrixToBinary(double **a, int m, int n, std::string filename) {
    std::ofstream file(filename, std::ios::binary|std::ios::out);
    for (int i = 0; i < m; i++) {
        file.write(reinterpret_cast<char *>(a[i]), n*sizeof(double));
    }
    file.close();
}

// write dimensions of n x m matrix to txt file
void writeMatrixDim(int n, int m, std::string filename) {
	std::ofstream file(filename + "_DIM.txt");
	file << n << " " << m;
	file.close();
}

template <typename myType>
myType **createAnyMatrix(int row, int col, myType value) {
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
myType *createAnyArray(myType value, int n) {
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

void destroyPointerMatrix(void **matr) {
	delete[](char *) matr[0];
	delete[] matr;

}

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum) {
	int j;
	long k;
	static long iy = 0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum = 1;
		else *idum = -(*idum);
		for (j = NTAB+7; j>=0; j--) {
			k = (*idum)/IQ;
			*idum = IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum)/IQ;
	*idum = IA*(*idum - k*IQ) - IR*k;
	if (*idum < 0) *idum += IM;
	j = iy/NDIV;
	iy = iv[j];
	iv[j] = *idum;
	if ((temp = (float)AM*iy) > RNMX) return (float)RNMX;
	else return temp;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX