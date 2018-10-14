#include <iostream>
#include <cmath>
#include <string>
#include <fstream>

double* createVector(double value, int n);

void integrateVerlet(int n, double beta, double FinalTime, double *pos_x, double *pos_y, double *vel_x, double *vel_y);

void integrateEuler(int n, double *pos_x, double *pos_y, double *vel_x, double *vel_y);

void doubleArrayToFile(double *v , int n, std::string filename);

double* linspace(double min, double max, int n);

void integrateVerletEJ2D(int n, double FinalTime, double *pos_Ex, double *pos_Ey, double *vel_Ex, double *vel_Ey,
double *pos_Jx, double *pos_Jy, double *vel_Jx, double *vel_Jy);

void integrateVerletEJ3D(int n, double FinalTime, double *pos_Ex, double *pos_Ey, double *pos_Ez, double *vel_Ex, double *vel_Ey,
double *vel_Ez, double *pos_Jx, double *pos_Jy, double *pos_Jz, double *vel_Jx, double *vel_Jy, double *vel_Jz);

void earthJupiterSun2D();

void earthJupiterSun3D();

void earthSun2D();

int main() {

    //earthJupiterSun3D();
    //earthJupiterSun2D();
    //earthSun2D();
    return 0;
}

void earthJupiterSun3D() {
    double *pos_Ex, *pos_Ey, *pos_Ez, *vel_Ex, *vel_Ey, *vel_Ez;
    double *pos_Jx, *pos_Jy, *pos_Jz, *vel_Jx, *vel_Jy, *vel_Jz;
    double FinalTime = 20.0; // years
    int n = 10000;

    pos_Ex = createVector(0.0, n);
    pos_Ey = createVector(0.0, n);
    pos_Ez = createVector(0.0, n);
    vel_Ex = createVector(0.0, n);
    vel_Ey = createVector(0.0, n);
    vel_Ez = createVector(0.0, n);

    pos_Jx = createVector(0.0, n);
    pos_Jy = createVector(0.0, n);
    pos_Jz = createVector(0.0, n);
    vel_Jx = createVector(0.0, n);
    vel_Jy = createVector(0.0, n);
    vel_Jz = createVector(0.0, n);

    pos_Ex[0] = 9.413801075750535e-01; // AU
    pos_Ey[0] = 3.379019986046322e-01; // AU
    pos_Ez[0] = -9.334104672733438e-05; // AU
    vel_Ex[0] = -5.994522787486753e-03*360.0; // AU/year
    vel_Ey[0] = 1.617377250092178e-02*360.0; // AU/year
    vel_Ez[0] = -1.732657683299539e-07*360; // AU/year

    pos_Jx[0] = -2.666952709077877; // AU
    pos_Jy[0] = -4.655671225645230; // AU
    pos_Jz[0] = 7.896515774211305e-02; // AU
    vel_Jx[0] = 6.458958874387921e-03*360.0; // AU/year
    vel_Jy[0] = -3.390642961368397e-03*360.0; // AU/year
    vel_Jz[0] = -1.303431975919576e-04*360.0; // AU/year


    integrateVerletEJ3D(n, FinalTime, pos_Ex, pos_Ey, pos_Ez, vel_Ex, vel_Ey, vel_Ez,
         pos_Jx, pos_Jy, pos_Jz, vel_Jx, vel_Jy, vel_Jz);

    doubleArrayToFile(pos_Ex, n, "pos_Ex");
    doubleArrayToFile(pos_Ey, n, "pos_Ey");
    doubleArrayToFile(pos_Ez, n, "pos_Ez");
    doubleArrayToFile(pos_Jx, n, "pos_Jx");
    doubleArrayToFile(pos_Jy, n, "pos_Jy");
    doubleArrayToFile(pos_Jz, n, "pos_Jz");

    delete[] pos_Ex;
    delete[] pos_Ey;
    delete[] pos_Ez;
    delete[] vel_Ex;
    delete[] vel_Ey;
    delete[] vel_Ez;
    delete[] pos_Jx;
    delete[] pos_Jy;
    delete[] pos_Jz;
    delete[] vel_Jx;
    delete[] vel_Jy;
    delete[] vel_Jz;
}

void earthJupiterSun2D() {
    double *pos_Ex, *pos_Ey, *vel_Ex, *vel_Ey;
    double *pos_Jx, *pos_Jy,*vel_Jx, *vel_Jy;
    double FinalTime = 20.0; // years
    int n = 10000;

    pos_Ex = createVector(0.0, n);
    pos_Ey = createVector(0.0, n);
    vel_Ex = createVector(0.0, n);
    vel_Ey = createVector(0.0, n);

    pos_Jx = createVector(0.0, n);
    pos_Jy = createVector(0.0, n);
    vel_Jx = createVector(0.0, n);
    vel_Jy = createVector(0.0, n);

    pos_Ex[0] = 9.413801075750535e-01; // AU
    pos_Ey[0] = 3.379019986046322e-01; // AU
    vel_Ex[0] = -5.994522787486753e-03*360.0; // AU/year
    vel_Ey[0] = 1.617377250092178e-02*360.0; // AU/year

    pos_Jx[0] = -2.666952709077877; // AU
    pos_Jy[0] = -4.655671225645230; // AU
    vel_Jx[0] = 6.458958874387921e-03*360.0; // AU/year
    vel_Jy[0] = -3.390642961368397e-03*360.0; // AU/year


    integrateVerletEJ2D(n, FinalTime, pos_Ex, pos_Ey, vel_Ex, vel_Ey,
         pos_Jx, pos_Jy, vel_Jx, vel_Jy);

    doubleArrayToFile(pos_Ex, n, "pos_Ex");
    doubleArrayToFile(pos_Ey, n, "pos_Ey");
    doubleArrayToFile(pos_Jx, n, "pos_Jx");
    doubleArrayToFile(pos_Jy, n, "pos_Jy");

    delete[] pos_Ex;
    delete[] pos_Ey;
    delete[] vel_Ex;
    delete[] vel_Ey;
    delete[] pos_Jx;
    delete[] pos_Jy;
    delete[] vel_Jx;
    delete[] vel_Jy;
}

void earthSun2D() {
    double *pos_x, *pos_y, *vel_x, *vel_y;
    double FinalTime = 8.0;
    int n = 100000;
    double beta = 3.0;
    //double ve = sqrt(8*acos(-1.0)*acos(-1.0)); // Escape velocity (analytic)

    pos_x = createVector(0.0, n);
    pos_y = createVector(0.0, n);
    vel_x = createVector(0.0, n);
    vel_y = createVector(0.0, n);

    pos_x[0] = 9.413801075750535e-01;
    pos_y[0] = 3.379019986046322e-01;
    vel_x[0] = -5.994522787486753e-03*360.0;
    vel_y[0] = 1.617377250092178E-02*360.0;


    integrateVerlet(n, beta, FinalTime, pos_x, pos_y, vel_x, vel_y);
    //integrateEuler(n, pos_x, pos_y, vel_x, vel_y);
    doubleArrayToFile(pos_x, n, "pos_x");
    doubleArrayToFile(pos_y, n, "pos_y");

    delete[] pos_x;
    delete[] pos_y;
    delete[] vel_x;
    delete[] vel_y;
}

void integrateVerlet(int n, double beta, double FinalTime, double *pos_x, double *pos_y, double *vel_x, double *vel_y) {
    double h = FinalTime/(n+1);
    double hh = h*h;
    double pi = 3.14159265359;
    double r, rn, a, an;
    double FourPi = 4*pi*pi;

    r = sqrt(pos_x[0]*pos_x[0] + pos_y[0]*pos_y[0]);
    a = -FourPi/pow(r, beta);

    for (int i = 0; i < (n-1); i++) {
        pos_x[i+1] = pos_x[i] + h*vel_x[i] + (hh/2.0)*a*pos_x[i];
        pos_y[i+1] = pos_y[i] + h*vel_y[i] + (hh/2.0)*a*pos_y[i];

        rn = sqrt(pos_x[i+1]*pos_x[i+1] + pos_y[i+1]*pos_y[i+1]);
        an = -FourPi/pow(rn, beta);

        vel_x[i+1] = vel_x[i] + (h/2.0)*(an*pos_x[i+1] + a*pos_x[i]);
        vel_y[i+1] = vel_y[i] + (h/2.0)*(an*pos_y[i+1] + a*pos_y[i]);

        r = rn;
        a = an;
    }
}

void integrateVerletEJ2D(int n, double FinalTime, double *pos_Ex, double *pos_Ey, double *vel_Ex, double *vel_Ey,
double *pos_Jx, double *pos_Jy, double *vel_Jx, double *vel_Jy) {
    double h = FinalTime/(n+1);
    double hh = h*h;
    double pi = 3.14159265359;
    double rE, rEn, rJ, rJn, rEJ, rEJn, aE, aJ, aEJ, aJE, aEn, aJn, aEJn, aJEn;
    double FourPi = 4*pi*pi;

    rE = sqrt(pos_Ex[0]*pos_Ex[0] + pos_Ey[0]*pos_Ey[0]);
    rJ = sqrt(pos_Jx[0]*pos_Jx[0] + pos_Jy[0]*pos_Jy[0]);
    rEJ = sqrt(pow((pos_Ex[0] - pos_Jx[0]), 2) + pow((pos_Ey[0] - pos_Jy[0]), 2));
    aE = -FourPi/pow(rE, 3);
    aJ = -FourPi/pow(rJ, 3);
    aEJ = -FourPi*0.00095/(pow(rEJ, 3));
    aJE = -FourPi*3e-3/(pow(rEJ, 3));


    for (int i = 0; i < (n-1); i++) {
        pos_Ex[i+1] = pos_Ex[i] + h*vel_Ex[i] + (hh/2.0)*(aE*pos_Ex[i] + aEJ*(pos_Ex[i] - pos_Jx[i]));
        pos_Ey[i+1] = pos_Ey[i] + h*vel_Ey[i] + (hh/2.0)*(aE*pos_Ey[i] + aEJ*(pos_Ey[i] - pos_Jy[i]));

        pos_Jx[i+1] = pos_Jx[i] + h*vel_Jx[i] + (hh/2.0)*(aJ*pos_Jx[i] + aJE*(pos_Jx[i] - pos_Ex[i]));
        pos_Jy[i+1] = pos_Jy[i] + h*vel_Jy[i] + (hh/2.0)*(aJ*pos_Jy[i] + aJE*(pos_Jy[i] - pos_Ey[i]));

        rEn = sqrt(pow(pos_Ex[i+1], 2.0) + pow(pos_Ey[i+1], 2.0));
        rJn = sqrt(pow(pos_Jx[i+1], 2.0) + pow(pos_Jy[i+1], 2.0));
        rEJn = sqrt(pow((pos_Ex[i+1] - pos_Jx[i+1]), 2) + pow((pos_Ey[i+1] - pos_Jy[i+1]), 2));
        aEn = -FourPi/pow(rEn, 3.0);
        aJn = -FourPi/pow(rJn, 3.0);
        aEJn = -FourPi*0.00095/(pow(rEJn, 3));
        aJEn = -FourPi*3e-3/(pow(rEJn, 3));

        vel_Ex[i+1] = vel_Ex[i] + (h/2.0)*((aEn*pos_Ex[i+1] + aEJn*(pos_Ex[i+1] - pos_Jx[i+1])) + (aE*pos_Ex[i] + aEJ*(pos_Ex[i] - pos_Jx[i])));
        vel_Ey[i+1] = vel_Ey[i] + (h/2.0)*((aEn*pos_Ey[i+1] + aEJn*(pos_Ey[i+1] - pos_Jy[i+1])) + (aE*pos_Ey[i] + aEJ*(pos_Ey[i] - pos_Jy[i])));

        vel_Jx[i+1] = vel_Jx[i] + (h/2.0)*((aJn*pos_Jx[i+1] + aJEn*(pos_Jx[i+1] - pos_Ex[i+1])) + (aJ*pos_Jx[i] + aJE*(pos_Jx[i] - pos_Ex[i])));
        vel_Jy[i+1] = vel_Jy[i] + (h/2.0)*((aJn*pos_Jy[i+1] + aJEn*(pos_Jy[i+1] - pos_Ey[i+1])) + (aJ*pos_Jy[i] + aJE*(pos_Jy[i] - pos_Ey[i])));

        rE = rEn;
        rJ = rJn;
        rEJ = rEJn;
        aE = aEn;
        aJ = aJn;
        aEJ = aEJn;
        aJE = aJEn;
    }

}

void integrateVerletEJ3D(int n, double FinalTime, double *pos_Ex, double *pos_Ey, double *pos_Ez, double *vel_Ex, double *vel_Ey,
double *vel_Ez, double *pos_Jx, double *pos_Jy, double *pos_Jz, double *vel_Jx, double *vel_Jy, double *vel_Jz) {
    double h = FinalTime/(n+1);
    double hh = h*h;
    double pi = 3.14159265359;
    double rE, rEn, rJ, rJn, rEJ, rEJn, aE, aJ, aEJ, aJE, aEn, aJn, aEJn, aJEn;
    double FourPi = 4*pi*pi;

    rE = sqrt(pow(pos_Ex[0], 2) + pow(pos_Ey[0],2) + pow(pos_Ez[0], 2));
    rJ = sqrt(pow(pos_Jx[0], 2) + pow(pos_Jy[0],2) + pow(pos_Jz[0], 2));
    rEJ = sqrt(pow((pos_Ex[0] - pos_Jx[0]), 2) + pow((pos_Ey[0] - pos_Jy[0]), 2) + pow((pos_Ez[0] - pos_Jz[0]), 2));
    aE = -FourPi/pow(rE, 3);
    aJ = -FourPi/pow(rJ, 3);
    aEJ = -FourPi*0.00095/(pow(rEJ, 3));
    aJE = -FourPi*3e-3/(pow(rEJ, 3));


    for (int i = 0; i < (n-1); i++) {
        pos_Ex[i+1] = pos_Ex[i] + h*vel_Ex[i] + (hh/2.0)*(aE*pos_Ex[i] + aEJ*(pos_Ex[i] - pos_Jx[i]));
        pos_Ey[i+1] = pos_Ey[i] + h*vel_Ey[i] + (hh/2.0)*(aE*pos_Ey[i] + aEJ*(pos_Ey[i] - pos_Jy[i]));
        pos_Ez[i+1] = pos_Ez[i] + h*vel_Ez[i] + (hh/2.0)*(aE*pos_Ez[i] + aEJ*(pos_Ez[i] - pos_Jz[i]));

        pos_Jx[i+1] = pos_Jx[i] + h*vel_Jx[i] + (hh/2.0)*(aJ*pos_Jx[i] + aJE*(pos_Jx[i] - pos_Ex[i]));
        pos_Jy[i+1] = pos_Jy[i] + h*vel_Jy[i] + (hh/2.0)*(aJ*pos_Jy[i] + aJE*(pos_Jy[i] - pos_Ey[i]));
        pos_Jz[i+1] = pos_Jz[i] + h*vel_Jz[i] + (hh/2.0)*(aJ*pos_Jz[i] + aJE*(pos_Jz[i] - pos_Ez[i]));

        rEn = sqrt(pow(pos_Ex[i+1], 2.0) + pow(pos_Ey[i+1], 2.0) + pow(pos_Ez[i+1], 2));
        rJn = sqrt(pow(pos_Jx[i+1], 2.0) + pow(pos_Jy[i+1], 2.0) + pow(pos_Jz[i+1], 2));
        rEJn = sqrt(pow((pos_Ex[i+1] - pos_Jx[i+1]), 2) + pow((pos_Ey[i+1] - pos_Jy[i+1]), 2) + pow((pos_Ez[i+1] - pos_Jz[i+1]), 2));
        aEn = -FourPi/pow(rEn, 3.0);
        aJn = -FourPi/pow(rJn, 3.0);
        aEJn = -FourPi*0.00095/(pow(rEJn, 3));
        aJEn = -FourPi*3e-3/(pow(rEJn, 3));

        vel_Ex[i+1] = vel_Ex[i] + (h/2.0)*((aEn*pos_Ex[i+1] + aEJn*(pos_Ex[i+1] - pos_Jx[i+1])) + (aE*pos_Ex[i] + aEJ*(pos_Ex[i] - pos_Jx[i])));
        vel_Ey[i+1] = vel_Ey[i] + (h/2.0)*((aEn*pos_Ey[i+1] + aEJn*(pos_Ey[i+1] - pos_Jy[i+1])) + (aE*pos_Ey[i] + aEJ*(pos_Ey[i] - pos_Jy[i])));
        vel_Ez[i+1] = vel_Ez[i] + (h/2.0)*((aEn*pos_Ez[i+1] + aEJn*(pos_Ez[i+1] - pos_Jz[i+1])) + (aE*pos_Ez[i] + aEJ*(pos_Ez[i] - pos_Jz[i])));

        vel_Jx[i+1] = vel_Jx[i] + (h/2.0)*((aJn*pos_Jx[i+1] + aJEn*(pos_Jx[i+1] - pos_Ex[i+1])) + (aJ*pos_Jx[i] + aJE*(pos_Jx[i] - pos_Ex[i])));
        vel_Jy[i+1] = vel_Jy[i] + (h/2.0)*((aJn*pos_Jy[i+1] + aJEn*(pos_Jy[i+1] - pos_Ey[i+1])) + (aJ*pos_Jy[i] + aJE*(pos_Jy[i] - pos_Ey[i])));
        vel_Jz[i+1] = vel_Jz[i] + (h/2.0)*((aJn*pos_Jz[i+1] + aJEn*(pos_Jz[i+1] - pos_Ez[i+1])) + (aJ*pos_Jz[i] + aJE*(pos_Jz[i] - pos_Ez[i])));

        rE = rEn;
        rJ = rJn;
        rEJ = rEJn;
        aE = aEn;
        aJ = aJn;
        aEJ = aEJn;
        aJE = aJEn;
    }

}

void integrateEuler(int n, double *pos_x, double *pos_y, double *vel_x, double *vel_y) {
    double h = 1.0/(n+1);
    double pi = 3.14159265359;
    double r;
    double FourPi = 4*pi*pi;

    for (int i = 0; i < (n-1); i++) {
        r = sqrt(pos_x[i]*pos_x[i] + pos_y[i]*pos_y[i]);

        pos_x[i+1] = pos_x[i] + h*vel_x[i];
        pos_y[i+1] = pos_y[i] + h*vel_y[i];

        vel_x[i+1] = vel_x[i] - h*FourPi/(r*r*r)*pos_x[i];
        vel_y[i+1] = vel_y[i] - h*FourPi /(r*r*r)*pos_y[i];
    }
}

double* createVector(double value, int n) {
	double *vector;
	vector = new double[n];
	for (int i = 0; i < n; i++) {
		vector[i] = value;
	}
	return vector;
}

void doubleArrayToFile(double *v , int n, std::string filename) {
	std::ofstream myfile(filename + ".txt");
	if (myfile.is_open()) {
		myfile << n << "\n";
		for (int i = 0; i < n; i++) {
			myfile << v[i] << "\n";
		}
	}
}

double* linspace(double min, double max, int n) {
    double *v;
    v = new double[n];
    double step = (max - min)/(n-1);
    v[0] = min;
	for (int i = 1; i < n; i++) {
        v[i] = min + i*step;
	}
    return v;
