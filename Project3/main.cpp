#include <string>
#include "NBS.h"
#include "planets_v2.h"
#include <iostream>

using namespace std;
using namespace StellarObjectsLibraryv2;

void earthSun();
void allPlanets();
void escapeVelocity();
void mercury();

int main() {

	mercury();
	//earthSun();
	//escapeVelocity();
	return 0;
}

void mercury() {
	int m, n;
	double finalTime;
	string filename;
	MassObject *planets;
	n = 1e9;

	m = 1;
	planets = new MassObject[m];
	MassObject RelMercury = { "RelMercury", Mercury.mass, 0.3075, 0, 0,
		0, 12.44, 0 };
	planets[0] = RelMercury;
	finalTime = 100 * 88.0 / 365.0;

	
	SolarSystemRelativistic mercury(planets, m, 0);
	mercury.setCenterMass(Sun.mass);
	mercury.verletSolveRel(finalTime, n, 0);
	mercury.perihelionPrecession(100);
	system("PAUSE");

}

void earthSun() {
	int m, n;
	int finalTime;
	string filename;
	MassObject *planets;

	m = 1;
	planets = new MassObject[m];

	MassObject Earth2D = { "2DEarth", Earth.mass, 1, 0, 0,
		0, 2.0*acos(-1.0), 0 };
	planets[0] = Earth2D;
	finalTime = 50;

	SolarSystem earth_sun(planets, m);
	earth_sun.setCenterMass(Sun.mass);
	filename = "earth_sun_euler_100";
	earth_sun.eulerSolve(finalTime, 100*finalTime);
	earth_sun.writeToFile(filename);
	filename = "earth_sun_verlet_100";
	earth_sun.verletSolve(finalTime, 100*finalTime);
	earth_sun.writeToFile(filename);
	earth_sun.setCenterMass(Sun.mass);
	filename = "earth_sun_euler_1000";
	earth_sun.eulerSolve(finalTime, 1000*finalTime);
	earth_sun.writeToFile(filename);
	filename = "earth_sun_verlet_1000";
	earth_sun.verletSolve(finalTime, 1000*finalTime);
	earth_sun.writeToFile(filename);
	earth_sun.setCenterMass(Sun.mass);
	filename = "earth_sun_euler_10000";
	earth_sun.eulerSolve(finalTime, 10000*finalTime);
	earth_sun.writeToFile(filename);
	filename = "earth_sun_verlet_10000";
	earth_sun.verletSolve(finalTime, 10000*finalTime);
	earth_sun.writeToFile(filename);
	filename = "earth_sun_euler_100000";
	earth_sun.eulerSolve(finalTime, 100000*finalTime);
	earth_sun.writeToFile(filename);
	filename = "earth_sun_verlet_100000";
	earth_sun.verletSolve(finalTime, 100000*finalTime);
	earth_sun.writeToFile(filename);
	delete[] planets;
}

void allPlanets() {
	int m, n;
	double finalTime;
	string filename;
	MassObject *planets;

	m = 10;
	planets = new MassObject[m];
	planets[0] = Sun;
	planets[1] = Mercury;
	planets[2] = Venus;
	planets[3] = Earth;
	planets[4] = Mars;
	planets[5] = Jupiter;
	planets[6] = Saturn;
	planets[7] = Uranus;
	planets[8] = Neptune;
	planets[9] = Pluto;
	n = 100000;
	finalTime = 250.0;
	filename = "whole";
	SolarSystem whole(planets, m);
	whole.verletSolve(finalTime, n);
	whole.writeToFile(filename);
	delete[] planets;
}

void escapeVelocity() {
    int m, n;
    double finalTime, *A;
    MassObject *planets;

    // part a) Earth-Sun system, fixed mass center.
    m = 1;
    planets = new MassObject[m];
    // what should the initial velocity be here?
    double pi = 3.14159265359;
    MassObject planet = {"Earth", 6e24, 1, 0, 0, 0, 0, 0};
    planets[0] = planet;
    finalTime = 3500.0;
    n = 1000000;


    double v1 = 8.875;
    planets[0].vy = v1;
    SolarSystem earth_sun(planets, m);
    earth_sun.setCenterMass(Sun.mass);
    earth_sun.verletSolve(finalTime, n);
    A = earth_sun.getAcceleration(1);
    double *R = earth_sun.getDistance(1);
    printf("a = %.20f, r = %.3f \n", A[0], R[0]);
    A = earth_sun.getAcceleration(n/2);
    R = earth_sun.getDistance(n/2);
    printf("a = %.20f, r = %.3f \n", A[0], R[0]);
    A = earth_sun.getAcceleration(n-1);
    R = earth_sun.getDistance(n-1);
    printf("a = %.20f, r = %.3f \n", A[0], R[0]);
    earth_sun.writeToFile("Earth_Sun_min");
    

    double v2 = 8.89;
    planets[0].vy = v2;
    earth_sun.verletSolve(finalTime, n);
    A = earth_sun.getAcceleration(1);
    R = earth_sun.getDistance(1);
    printf("\na = %.20f, r = %.3f \n", A[0], R[0]);
    A = earth_sun.getAcceleration(n/2);
    R = earth_sun.getDistance(n/2);
    printf("a = %.20f, r = %.3f \n", A[0], R[0]);
    A = earth_sun.getAcceleration(n-1);
    R = earth_sun.getDistance(n-1);
    printf("a = %.20f, r = %.3f \n", A[0], R[0]);
    earth_sun.writeToFile("Earth_Sun_max");

    double ve = (v1 + v2)/2;
    double ana = sqrt(8)*pi;
    printf("\nAnalytical %.20f \n", ana);
    printf("Ve %.20f \n", ve);
    printf("relative error %.5f %% \n", fabs(ana - ve)/ana*100.0);

    finalTime = 20.0;
    n = 10000;

    v1 = 6.8;
    v2 = 7.4;
    double v3 = 8.0;
    double v4 = 8.8;
    double beta;
    string filename;
    for (int i = 0; i < 4; i++) {
        beta = 2 + (double) i/3.0;
        printf("\nBeta = %.3f\n", beta);
        filename = "Beta_" + to_string(i) + "_";
        earth_sun.setBeta(beta);

        printf("v = %.3f\n", v1);
        planets[0].vy = v1;
        earth_sun.verletSolve(finalTime, n);
        A = earth_sun.getAcceleration(n-1);
        R = earth_sun.getDistance(n-1);
        printf("a = %.20f, r = %.3f \n", A[0], R[0]);
        earth_sun.writeToFile(filename + to_string(1));

        printf("v = %.3f\n", v2);
        planets[0].vy = v2;
        earth_sun.verletSolve(finalTime, n);
        A = earth_sun.getAcceleration(n-1);
        R = earth_sun.getDistance(n-1);
        printf("a = %.20f, r = %.3f \n", A[0], R[0]);
        earth_sun.writeToFile(filename + to_string(2));

        printf("v = %.3f\n", v3);
        planets[0].vy = v3;
        earth_sun.verletSolve(finalTime, n);
        A = earth_sun.getAcceleration(n-1);
        R = earth_sun.getDistance(n-1);
        printf("a = %.20f, r = %.3f \n", A[0], R[0]);
        earth_sun.writeToFile(filename + to_string(3));


        printf("v = %.3f\n", v4);
        planets[0].vy = v4;
        earth_sun.verletSolve(finalTime, n);
        A = earth_sun.getAcceleration(n-1);
        R = earth_sun.getDistance(n-1);
        printf("a = %.20f, r = %.3f \n", A[0], R[0]);
        earth_sun.writeToFile(filename + to_string(4));
    }

    delete[] planets;
    delete[] R;
    delete[] A;
}
