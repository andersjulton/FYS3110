#include <iostream>
#include <cmath>

double Z(double temp) {
	double Z;
	Z = 2*exp(8/temp) + 2*exp(-8/temp) + 12;
	return Z;
}

double eMeanA(double temp) {
	double E;
	E = -(1.0/Z(temp))*(16*exp(8/temp) - 16*exp(-8/temp));
	return E;
}

double eSqMeanA(double temp) {
	double ESq;
	ESq = (1.0/Z(temp))*(128*exp(8/temp) + 128*exp(-8/temp));
	return ESq;
}

double mabsMeanA(double temp) {
	double M;
	M = (1.0/Z(temp))*(8*exp(8/temp) + 16);
	return M;
}

double mSqMeanA(double temp) {
	double MSq;
	MSq = (1.0/Z(temp))*(32*exp(8/temp) + 32);
	return MSq;
}

double CVA(double temp) {
	double CV;
	CV = (1.0/(temp*temp))*(eSqMeanA(temp) - pow(eMeanA(1.0), 2));
	return CV;
}

double suscA(double temp) {
	double susc;
	susc = (1.0/(temp))*(mSqMeanA(temp) - mabsMeanA(temp)*mabsMeanA(temp));
	return susc;
}

void expValues(double temp) {
	std::cout << "\nAnalytical values for temperature = " << temp << '\n';
	std::cout << "Mean energy = " << eMeanA(temp) << '\n';
	std::cout << "Mean energy squared = " << eSqMeanA(temp) << '\n';
	std::cout << "Mean absolute magnetization = " << mabsMeanA(temp) << '\n';
	std::cout << "Mean magnetization squared = " << mSqMeanA(temp) << '\n';
	std::cout << "Specific heat = " << CVA(temp) << '\n';
	std::cout << "Magnetic susceptibility = " << suscA(temp) << '\n';
}