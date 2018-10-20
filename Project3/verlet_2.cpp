	#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <stdlib.h>

struct MassObject {
	std::string name;
	double mass;
	// position:
	double x;
	double y;
	double z;
	//velocity:
	double vx;
	double vy;
	double vz;
};

namespace StellarObjectsLibrary {

    struct MassObject Sun = {"Sun", 2e30, -1.474798697566110E-04, 7.248922398410132E-03, -7.268118257496195E-05,
        -7.585746627696994E-06*365.0, 2.590655894257502E-06*365.0, 1.896113821987181E-07*365.0};

    struct MassObject Mercury = {"Mercury", 3.3e23, -1.534743586808411E-01, -4.321686982270908E-01, -2.191308023125118E-02,
         2.090296530646647E-02*365.0, -7.871323434203551E-03*365.0, -2.561512430419179E-03*365.0};

    struct MassObject Venus = {"Venus", 4.9e24, 7.066148066072473E-01, 1.672271996665292E-01, -3.866246273128392E-02,
        -4.546772256640751E-03*365.0, 1.963882295765339E-02*365.0, 5.315229360444450E-04*365.0};

    struct MassObject Earth = {"Earth", 6e24, 9.413801075750535e-01, 3.379019986046322E-01, -9.334104672733438E-05,
        -5.994522787486753E-03*365.0, 1.617377250092178E-02*365.0, -1.732657683299539E-07*365.0};

    struct MassObject Mars = {"Mars", 6.6e23, 1.377524608498332, -1.476563536087052E-01, -3.712365888122099E-02,
         2.090430895774286E-03*365.0, 1.510431174964668E-02*365.0, 2.651531936464447E-04*365.0};

    struct MassObject Jupiter = {"Jupiter", 1.9e27, -2.666952709077877, -4.655671225645230, 7.896515774211305E-02,
        6.458958874387921E-03*365.0, -3.390642961368397E-03*365.0, -1.303431975919576E-04*365.0};

    struct MassObject Saturn = {"Saturn", 5.5e26, 1.549159416633817, -9.935197133379326, 1.110808051816651E-01,
         5.204592319579984E-03*365.0, 8.422049097583516E-04*365.0, -2.220199495831516E-04*365.0};

    struct MassObject Uranus = {"Uranus", 8.8e25, 1.717591494590517E+01, 9.997664143031971, -1.853846122526302E-01,
         -2.007356242188686E-03*365.0, 3.215850240122884E-03*365.0, 3.800786256690271E-05*365.0};

    struct MassObject Neptune = {"Neptune", 1.03e26, 2.892029941220658E+01, -7.722539840090450, -5.074661608355547E-01,
         7.890393808745537E-04*365.0, 3.051931545808817E-03*365.0, -8.068388455453793E-05*365.0};

    struct MassObject Pluto = {"Pluto", 1.31e22, 1.164390186496279E+01, -3.157511878129099E+01, 1.062859645894982E-02,
         3.018604420452015E-03*365.0, 4.214379702145380E-04*365.0, -9.301537706126110E-04*365.0};
}

double* createVector(double value, int n);

double* linspace(double min, double max, int n);

double **createMatrix(int m, int n);

void deleteMatrix(double **mat, int n);

void doubleMatrixToFile(double **v , int n, int m, std::string filename);

void doubleArrayToFile(double *v , int n, std::string filename);

void integrateEuler(int n, double *pos_x, double *pos_y, double *vel_x, double *vel_y);

void integrateVerlet(int n, double beta, double mass, double finalTime, double *pos_x, double *pos_y, double *pos_z,
double *vel_x, double *vel_y, double *vel_z, bool testEnergy);

void integrateVerletEJ3D(int n, double finalTime, double *pos_Ex, double *pos_Ey, double *pos_Ez, double *vel_Ex, double *vel_Ey,
double *vel_Ez, double *pos_Jx, double *pos_Jy, double *pos_Jz, double *vel_Jx, double *vel_Jy, double *vel_Jz);

void integrateVerletPlanetsStruct(int n, double **pos_x, double **pos_y, double **pos_z, double **vel_x, double **vel_y, double **vel_z,
double finalTime, struct MassObject *planets, int m);

void mercuryIntegrate(int n, bool rel, double finalTime, double *pos_x, double *pos_y, double *vel_x, double *vel_y);

void earthSun();

void earthJupiterSun3D();

void planets();

void mercury(bool rel);

void planetsStruct();

int main() {

	mercury(false);
	mercury(true);
    return 0;
}

void earthSun() {
    int n = 10000;
    double finalTime = 20.0;
	double pi = 3.14159265359;
    double *pos_x, *pos_y, *pos_z, *vel_x, *vel_y, *vel_z;
    pos_x = createVector(0, n);
    pos_y = createVector(0, n);
    pos_z = createVector(0, n);
    vel_x = createVector(0, n);
    vel_y = createVector(0, n);
    vel_z = createVector(0, n);

	pos_x[0] = 1.0;
	pos_y[0] = 0.0;
	pos_z[0] = 0.0;
	vel_x[0] = 0.0;
	vel_y[0] = sqrt(8*pi*pi);
	vel_z[0] = 0.0;

    integrateVerlet(n, 3, StellarObjectsLibrary::Earth.mass, finalTime,
		pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, false);

    doubleArrayToFile(pos_x, n, "ES_pos_x");
    doubleArrayToFile(pos_y, n, "ES_pos_y");
    doubleArrayToFile(pos_z, n, "ES_pos_z");

    delete[] pos_x;
    delete[] pos_y;
    delete[] pos_z;
    delete[] vel_x;
    delete[] vel_y;
    delete[] vel_z;
}

void mercury(bool rel) {
	int n = 1e7;
    double finalTime = 100.0*365.0/88.0;
	double pi = 3.14159265359;
    double *pos_x, *pos_y, *vel_x, *vel_y;
    pos_x = createVector(0, n);
    pos_y = createVector(0, n);
    vel_x = createVector(0, n);
    vel_y = createVector(0, n);

	pos_x[0] = 0.3075;
	pos_y[0] = 0.0;
	vel_x[0] = 0.0;
	vel_y[0] = 12.44;

	int N = n - n/100;
	double r;
	double rmin;
	int index;

	mercuryIntegrate(n, rel, finalTime, pos_x, pos_y, vel_x, vel_y);

	rmin = sqrt(pow(pos_x[N-1], 2) + pow(pos_y[N-1], 2));
	for (int i = N; i < n; i++) {
		r = sqrt(pow(pos_x[i], 2) + pow(pos_y[i], 2));
		if (r < rmin) {
			rmin = r;
			index = i;
		}
	}

	std::cout << "Position of perihelion is: (" << pos_x[index] << "," << pos_y[index] << ")\n";
	std::cout << "Argument of perihelion is: " << atan2(pos_y[index], pos_x[index])*206264.806 << " arcseconds/century.\n";

    delete[] pos_x;
    delete[] pos_y;
    delete[] vel_x;
    delete[] vel_y;
}

void planetsStruct() {
    int n = 100000;
    double finalTime = 250.0;
    int m = 10;
    double **pos_x, **pos_y, **pos_z, **vel_x, **vel_y, **vel_z;
    pos_x = createMatrix(m, n);
    pos_y = createMatrix(m, n);
    pos_z = createMatrix(m, n);
    vel_x = createMatrix(m, n);
    vel_y = createMatrix(m, n);
    vel_z = createMatrix(m, n);

    struct MassObject *planets;
    planets = new struct MassObject[m];

    planets[0] = StellarObjectsLibrary::Sun;
    planets[1] = StellarObjectsLibrary::Mercury;
    planets[2] = StellarObjectsLibrary::Venus;
    planets[3] = StellarObjectsLibrary::Earth;
    planets[4] = StellarObjectsLibrary::Mars;
    planets[5] = StellarObjectsLibrary::Jupiter;
    planets[6] = StellarObjectsLibrary::Saturn;
    planets[7] = StellarObjectsLibrary::Uranus;
    planets[8] = StellarObjectsLibrary::Neptune;
    planets[9] = StellarObjectsLibrary::Pluto;

    integrateVerletPlanetsStruct(n, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, finalTime, planets, m);

    doubleMatrixToFile(pos_x, n, m, "pos_x");
    doubleMatrixToFile(pos_y, n, m, "pos_y");
    doubleMatrixToFile(pos_z, n, m, "pos_z");

    deleteMatrix(pos_x, n);
    deleteMatrix(pos_y, n);
    deleteMatrix(pos_z, n);
    deleteMatrix(vel_x, n);
    deleteMatrix(vel_y, n);
    deleteMatrix(vel_z, n);
}

void earthJupiterSun3D() {
    double *pos_Ex, *pos_Ey, *pos_Ez, *vel_Ex, *vel_Ey, *vel_Ez;
    double *pos_Jx, *pos_Jy, *pos_Jz, *vel_Jx, *vel_Jy, *vel_Jz;
    double finalTime = 20.0; // years
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
    vel_Ex[0] = -5.994522787486753e-03*365.0; // AU/year
    vel_Ey[0] = 1.617377250092178e-02*365.0; // AU/year
    vel_Ez[0] = -1.732657683299539e-07*360; // AU/year

    pos_Jx[0] = -2.666952709077877; // AU
    pos_Jy[0] = -4.655671225645230; // AU
    pos_Jz[0] = 7.896515774211305e-02; // AU
    vel_Jx[0] = 6.458958874387921e-03*365.0; // AU/year
    vel_Jy[0] = -3.390642961368397e-03*365.0; // AU/year
    vel_Jz[0] = -1.303431975919576e-04*365.0; // AU/year


    integrateVerletEJ3D(n, finalTime, pos_Ex, pos_Ey, pos_Ez, vel_Ex, vel_Ey, vel_Ez,
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

void mercuryIntegrate(int n, bool rel, double finalTime, double *pos_x, double *pos_y, double *vel_x, double *vel_y) {
    double h = finalTime/n;
    double hh = h*h;
    double pi = 3.14159265359;
	double c = 63239.7263; // Speed of light AU/yr
	double cc = c*c;
    double r, rn, a, an, l, ll;
    double fourPiPi = 4*pi*pi;

	r = sqrt(pow(pos_x[0], 2) + pow(pos_y[0], 2));

	if (rel) {l = pos_x[0]*vel_y[0] - pos_y[0]*vel_x[0];}
	else {l = 0;}

	ll = l*l;
    a = -fourPiPi*(1.0 + 3.0*ll/(r*r*cc))/pow(r, 3);
	std::cout << "Angular momentum before simulation : " << l << '\n';

    for (int i = 0; i < (n-1); i++) {
        pos_x[i+1] = pos_x[i] + h*vel_x[i] + (hh/2.0)*a*pos_x[i];
		pos_y[i+1] = pos_y[i] + h*vel_y[i] + (hh/2.0)*a*pos_y[i];

        rn = sqrt(pos_x[i+1]*pos_x[i+1] + pos_y[i+1]*pos_y[i+1]);
        an = -fourPiPi*(1.0 + 3.0*ll/(rn*rn*cc))/pow(rn, 3);

        vel_x[i+1] = vel_x[i] + (h/2.0)*(an*pos_x[i+1] + a*pos_x[i]);
        vel_y[i+1] = vel_y[i] + (h/2.0)*(an*pos_y[i+1] + a*pos_y[i]);

        r = rn;
        a = an;
    }
	std::cout << "Angular momentum after simulation : " << pos_x[n-1]*vel_y[n-1] - pos_y[n-1]*vel_x[n-1] << '\n';
}

void integrateVerlet(int n, double beta, double mass, double finalTime, double *pos_x, double *pos_y, double *pos_z,
	double *vel_x, double *vel_y, double *vel_z, bool testEnergy) {
    double h = finalTime/(n+1);
	double count = 0;
	double U0, Un, L0, Ln;
    double hh = h*h;
    double pi = 3.14159265359;
    double r, rn, a, an, v;
    double fourPiPi = 4*pi*pi;

	r = sqrt(pow(pos_x[0], 2) + pow(pos_y[0], 2) + pow(pos_z[0], 2));
    v = sqrt(pow(vel_x[0], 2) + pow(vel_y[0], 2) + pow(vel_z[0], 2));
    a = -fourPiPi/pow(r, beta);
	U0 = -fourPiPi/r;
	L0 = r*v;

    for (int i = 0; i < (n-1); i++) {
        pos_x[i+1] = pos_x[i] + h*vel_x[i] + (hh/2.0)*a*pos_x[i];
		pos_y[i+1] = pos_y[i] + h*vel_y[i] + (hh/2.0)*a*pos_y[i];
        pos_z[i+1] = pos_z[i] + h*vel_z[i] + (hh/2.0)*a*pos_z[i];

        rn = sqrt(pos_x[i+1]*pos_x[i+1] + pos_y[i+1]*pos_y[i+1] + pos_z[i+1]*pos_z[i+1]);
        an = -fourPiPi/pow(rn, beta);

        vel_x[i+1] = vel_x[i] + (h/2.0)*(an*pos_x[i+1] + a*pos_x[i]);
        vel_y[i+1] = vel_y[i] + (h/2.0)*(an*pos_y[i+1] + a*pos_y[i]);
		vel_z[i+1] = vel_z[i] + (h/2.0)*(an*pos_z[i+1] + a*pos_z[i]);
		v = sqrt(pow(vel_x[i+1], 2) + pow(vel_y[i+1], 2) + pow(vel_z[i+1], 2));

        r = rn;
        a = an;
		count++;
		if (0.5*v*v > fourPiPi/r) {
			std::cout << "Earth escaped the Sun with speed = " << v << " AU/yr\n";
			exit (EXIT_FAILURE);
		}
		if (testEnergy) {
			if (count = n/100.0) {
				Un = -fourPiPi/r;
				Ln = r*v;
				if (fabs(Un - U0) > 1e-2) {
					std::cout << "Potential energy is not conserved" << '\n';
					std::cout << "Un = " << Un << " U0 = " << U0 << '\n';
					exit (EXIT_FAILURE);
				}
				if (fabs(Ln - L0) > 1e-2) {
					std::cout << "Angluar momentum is not conserved" << '\n';
					std::cout << "Ln = " << Ln << " L0 = " << L0 << '\n';
					exit (EXIT_FAILURE);
				}
				count = 0;
			}
		}
    }
}

void integrateVerletPlanetsStruct(int n, double **pos_x, double **pos_y, double **pos_z, double **vel_x, double **vel_y, double **vel_z,
double finalTime, MassObject *planets, int m) { // NEW FORMAT WITH STRUCT
    double h = finalTime/(n+1);
    double hh = h*h;
    double pi = 3.14159265359;
    double fourPiPi = 4*pi*pi;
    double ax, ay, az, r, denum;
    double SM = 2e30;
    double **A;

    A = createMatrix(m, 3);
    for (int i = 0; i < m; i++) {
        ax = 0;
        ay = 0;
        az = 0;
        for (int j = 0; j < m; j++) {
            r = sqrt(pow(planets[i].x - planets[j].x, 2.0)
            + pow(planets[i].y - planets[j].y, 2.0) + pow(planets[i].z - planets[j].z, 2.0));
            if (r == 0) {
                continue;
            }
            denum = SM*pow(r, 3);
            ax -= fourPiPi*(planets[i].x - planets[j].x)*planets[j].mass/denum;
            ay -= fourPiPi*(planets[i].y - planets[j].y)*planets[j].mass/denum;
            az -= fourPiPi*(planets[i].z - planets[j].z)*planets[j].mass/denum;
        }
        pos_x[i][0] = planets[i].x;
        pos_y[i][0] = planets[i].y;
        pos_z[i][0] = planets[i].z;
        vel_x[i][0] = planets[i].vx;
        vel_y[i][0] = planets[i].vy;
        vel_z[i][0] = planets[i].vz;
        A[i][0] = ax;
        A[i][1] = ay;
        A[i][2] = az;

    }
    for (int i = 0; i < n-1; i++) {
        for (int j = 0; j < m; j++) {
            pos_x[j][i+1] = pos_x[j][i] + h*vel_x[j][i] + (hh/2.0)*A[j][0];
            pos_y[j][i+1] = pos_y[j][i] + h*vel_y[j][i] + (hh/2.0)*A[j][1];
            pos_z[j][i+1] = pos_z[j][i] + h*vel_z[j][i] + (hh/2.0)*A[j][2];
        }
        for (int j = 0; j < m; j++) {
            ax = 0;
            ay = 0;
            az = 0;
            for (int k = 0; k < m; k++) {
                r = sqrt(pow(pos_x[j][i+1] - pos_x[k][i+1], 2.0)  + pow(pos_y[j][i+1] - pos_y[k][i+1], 2.0) + pow(pos_z[j][i+1] - pos_z[k][i+1], 2.0));
                if (r == 0) {
                    continue;
                }
                denum = SM*pow(r, 3);
                ax -= fourPiPi*(pos_x[j][i+1] - pos_x[k][i+1])*planets[k].mass/denum;
                ay -= fourPiPi*(pos_y[j][i+1] - pos_y[k][i+1])*planets[k].mass/denum;
                az -= fourPiPi*(pos_z[j][i+1] - pos_z[k][i+1])*planets[k].mass/denum;
            }
            vel_x[j][i+1] = vel_x[j][i] + (h/2.0)*(A[j][0] + ax);
            vel_y[j][i+1] = vel_y[j][i] + (h/2.0)*(A[j][1] + ay);
            vel_z[j][i+1] = vel_z[j][i] + (h/2.0)*(A[j][2] + az);

            A[j][0] = ax;
            A[j][1] = ay;
            A[j][2] = az;
        }
    }
    deleteMatrix(A, 3);
}

void integrateVerletEJ3D(int n, double finalTime, double *pos_Ex, double *pos_Ey, double *pos_Ez, double *vel_Ex, double *vel_Ey,
double *vel_Ez, double *pos_Jx, double *pos_Jy, double *pos_Jz, double *vel_Jx, double *vel_Jy, double *vel_Jz) {
    double h = finalTime/(n+1);
    double hh = h*h;
    double pi = 3.14159265359;
    double rE, rEn, rJ, rJn, rEJ, rEJn, aE, aJ, aEJ, aJE, aEn, aJn, aEJn, aJEn;
    double fourPiPi = 4*pi*pi;

    rE = sqrt(pow(pos_Ex[0], 2) + pow(pos_Ey[0], 2) + pow(pos_Ez[0], 2));
    rJ = sqrt(pow(pos_Jx[0], 2) + pow(pos_Jy[0], 2) + pow(pos_Jz[0], 2));
    rEJ = sqrt(pow((pos_Ex[0] - pos_Jx[0]), 2) + pow((pos_Ey[0] - pos_Jy[0]), 2) + pow((pos_Ez[0] - pos_Jz[0]), 2));
    aE = -fourPiPi/pow(rE, 3);
    aJ = -fourPiPi/pow(rJ, 3);
    aEJ = -fourPiPi*0.00095/(pow(rEJ, 3));
    aJE = -fourPiPi*3e-3/(pow(rEJ, 3));

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
        aEn = -fourPiPi/pow(rEn, 3.0);
        aJn = -fourPiPi/pow(rJn, 3.0);
        aEJn = -fourPiPi*0.00095/(pow(rEJn, 3));
        aJEn = -fourPiPi*3e-3/(pow(rEJn, 3));

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
    double fourPiPi = 4*pi*pi;

    for (int i = 0; i < (n-1); i++) {
        r = sqrt(pos_x[i]*pos_x[i] + pos_y[i]*pos_y[i]);

        pos_x[i+1] = pos_x[i] + h*vel_x[i];
        pos_y[i+1] = pos_y[i] + h*vel_y[i];

        vel_x[i+1] = vel_x[i] - h*fourPiPi/(r*r*r)*pos_x[i];
        vel_y[i+1] = vel_y[i] - h*fourPiPi /(r*r*r)*pos_y[i];
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

void doubleMatrixToFile(double **v , int n, int m, std::string filename) {
	std::ofstream myfile(filename + ".txt");
	if (myfile.is_open()) {
		myfile << n << " " <<  m << "\n";
		for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                myfile << v[j][i] << " ";
            }
            myfile << '\n';
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
}

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

void deleteMatrix(double **mat, int n) {
	for (int i = 0; i < n; i++) {
		delete[] mat[i];
	}
	delete[] mat;
}
