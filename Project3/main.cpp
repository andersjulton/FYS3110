#include <string>
#include "NBS.h"
#include <cmath>
#include "planetsLib.h"
#include <iostream>

using namespace std;
using namespace StellarObjectsLibrary;

void earthSun();
void allPlanets();
void escapeVelocity();
void mercury();
void earthJupiter();
void sunEarthJupiter();
void time();

int main() {

	//sunEarthJupiter();
	//earthJupiter();
	//earthJupiter_mass();
	//mercury();
	earthSun();
	//allPlanets();
	//escapeVelocity();
	//time();
	
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
	finalTime = 100;

	SolarSystemRelativistic mercury(planets, m);
	mercury.setCenterMass(Sun.mass);
	mercury.perihelionPrecession(0, finalTime, n, 100);
	system("PAUSE");
	mercury.destroy();
	delete[] planets;
}

void time() {
	int m = 1;
	MassObject *planets = new MassObject[m];
	MassObject Earth2D = { "2DEarth", Earth.mass, 1, 0, 0, 0, 2.0*acos(-1.0), 0 };
	planets[0] = Earth2D;

    SolarSystem earth_sun(planets, m);
    earth_sun.setCenterMass(1.0);

	double finalTime = 50;
    int n = 1000;
    double timeE, timeV;
    for (int i = 0; i < 4; i++) {
    	printf("n = %d\n", n);
        timeE = earth_sun.timeEulerSolve(finalTime, (int) (n*finalTime));
        printf("euler = %f\n", timeE);
        timeV = earth_sun.timeVerletSolve(finalTime, (int) (n*finalTime));
        printf("verlet = %f\n", timeV);
        printf("ratio = %f\n", timeV/timeE);
        n = n*10;
    }
    earth_sun.destroy();
	delete[] planets;

}

void earthSun() {
	string filename1 = "earth_sun_euler_";
	string filename2 = "earth_sun_verlet_";
	MassObject *planets;

	int m = 1;
	int points = 100;
	planets = new MassObject[m];

	MassObject Earth2D = { "2DEarth", Earth.mass, 1, 0, 0, 0, 2.0*acos(-1.0), 0 };
	planets[0] = Earth2D;

    SolarSystem earth_sun(planets, m);
    earth_sun.setCenterMass(1.0);

	double finalTime = 50;
    int n = 100;
    for (int i = 0; i < 4; i++) {
        //earth_sun.eulerSolve(finalTime, (int) (n*finalTime));
        //earth_sun.writeToFile(filename1 + to_string(n));
		//earth_sun.conservation(filename1 + to_string(n), points);

        //earth_sun.verletSolve(finalTime, (int) (n*finalTime));
        //earth_sun.writeToFile(filename2 + to_string(n));
		//earth_sun.conservation(filename2 + to_string(n), points);
        n = n*10;
    }
	earth_sun.setBeta(3.0);
	earth_sun.verletSolve(finalTime, 1e6);
	earth_sun.writeToFile(filename2 + to_string(n));
    earth_sun.destroy();
	delete[] planets;
}

void sunEarthJupiter() {
	int m, n;
	double finalTime = 50.0;
	double Rx, Ry, Rz, M;
	Rx = 0;
	Ry = 0;
	Rz = 0;
	M = 0;
	string filename = "sunEarthJupiter";
	MassObject *planets;

	n = 1e5;
	m = 3;

	planets = new MassObject[m];
	planets[0] = Sun;
	planets[1] = Earth;
	planets[2] = Jupiter;
	
	for (int i = 0; i < m; i++) {
		M += planets[i].mass;
		Rx += planets[i].mass*planets[i].x;
		Ry += planets[i].mass*planets[i].y;
		Rz += planets[i].mass*planets[i].z;
	}
	for (int i = 0; i < m; i++) {
		planets[i].x -= Rx/M;
		planets[i].y -= Ry/M;
		planets[i].z -= Rz/M;
	}
	planets[0].vx = -(planets[1].vx*planets[1].mass + planets[2].vx*planets[2].mass);
	planets[0].vy = -(planets[1].vy*planets[1].mass + planets[2].vy*planets[2].mass);
	planets[0].vz = -(planets[1].vz*planets[1].mass + planets[2].vz*planets[2].mass);
	
	SolarSystem sunEarthJupiter(planets, m);
	sunEarthJupiter.verletSolve(finalTime, n);
	sunEarthJupiter.writeToFile(filename);

	sunEarthJupiter.destroy();
	delete[] planets;
}

void earthJupiter() {
	int m, n, points;
	double finalTime;
	string filename;
	MassObject *planets;

	n = 1e6;
	m = 2;
	finalTime = 50.0;
	points = 1000;
	filename = "earthJupiter";

	planets = new MassObject[m];
	MassObject earth, jupiter;
	earth.mass = Earth.mass;
	jupiter.mass = Jupiter.mass;

	// Positions with sun at center
	earth.x = 9.415200029562484E-01;
	earth.y = 3.306556713288177E-01;
	earth.z = -2.047028652191917E-05;
	earth.vx = -5.986939558621437E-03*365.0;
	earth.vy = 1.617117290454856E-02*365.0;
	earth.vz = -3.628095647481011E-07*365.0;
	
	jupiter.x = -2.666812813696682;
	jupiter.y = -4.662917552921043;
	jupiter.z = 7.903802850231831E-02;
	jupiter.vx = 6.466542103253237E-03*365.0;
	jupiter.vy = -3.393242557741619E-03*365.0;
	jupiter.vz = -1.305327413883757E-04*365.0;

	planets[0] = earth;
	planets[1] = jupiter;

	SolarSystem earthSun(planets, 1);
	earthSun.setCenterMass(Sun.mass);
	earthSun.verletSolve(finalTime, n);
	//earthSun.writeToFile("earthSun");
	earthSun.conservation("earthSun", points);

	SolarSystem earthJupiter(planets, m);
	earthJupiter.setCenterMass(Sun.mass);
	earthJupiter.verletSolve(finalTime, n);
	//earthJupiter.writeToFile(filename);
	earthJupiter.conservation(filename, points);

	planets[1].mass = Jupiter.mass*10.0;
	SolarSystem earthJupiter10(planets, m);
	earthJupiter10.setCenterMass(Sun.mass);
	earthJupiter10.verletSolve(finalTime, n);
	//earthJupiter10.writeToFile(filename + "10x");
	earthJupiter10.conservation(filename+"10x", points);

	planets[1].mass = Jupiter.mass*1000.0;
	SolarSystem earthJupiter1000(planets, m);
	earthJupiter1000.setCenterMass(Sun.mass);
	earthJupiter1000.verletSolve(finalTime, n);
	//earthJupiter1000.writeToFile(filename + "1000x");
	earthJupiter1000.conservation(filename+"1000x", points);

	earthSun.destroy();
	earthJupiter.destroy();
	earthJupiter10.destroy();
	earthJupiter1000.destroy();
	
	delete[] planets;
}

void allPlanets() {
	int m, n, points;
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

	n = 1e6;
	finalTime = 1000.0;
	points = 1000;
	filename = "whole";
	SolarSystem whole(planets, m);
	whole.verletSolve(finalTime, n);
	whole.writeToFile(filename);
	whole.conservation(filename, points);

	
    whole.destroy();
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
    delete[] R; delete[] A;
    A = earth_sun.getAcceleration(n/2);
    R = earth_sun.getDistance(n/2);
    printf("a = %.20f, r = %.3f \n", A[0], R[0]);
    delete[] R; delete[] A;
    A = earth_sun.getAcceleration(n-1);
    R = earth_sun.getDistance(n-1);
    printf("a = %.20f, r = %.3f \n", A[0], R[0]);
    delete[] R; delete[] A;
    earth_sun.writeToFile("Earth_Sun_min");


    double v2 = 8.89;
    planets[0].vy = v2;
    earth_sun.verletSolve(finalTime, n);
    A = earth_sun.getAcceleration(1);
    R = earth_sun.getDistance(1);
    printf("\na = %.20f, r = %.3f \n", A[0], R[0]);
    delete[] R; delete[] A;
    A = earth_sun.getAcceleration(n/2);
    R = earth_sun.getDistance(n/2);
    printf("a = %.20f, r = %.3f \n", A[0], R[0]);
    delete[] R; delete[] A;
    A = earth_sun.getAcceleration(n-1);
    R = earth_sun.getDistance(n-1);
    printf("a = %.20f, r = %.3f \n", A[0], R[0]);
    delete[] R; delete[] A;
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
        delete[] R; delete[] A;

        printf("v = %.3f\n", v2);
        planets[0].vy = v2;
        earth_sun.verletSolve(finalTime, n);
        A = earth_sun.getAcceleration(n-1);
        R = earth_sun.getDistance(n-1);
        printf("a = %.20f, r = %.3f \n", A[0], R[0]);
        earth_sun.writeToFile(filename + to_string(2));
        delete[] R; delete[] A;

        printf("v = %.3f\n", v3);
        planets[0].vy = v3;
        earth_sun.verletSolve(finalTime, n);
        A = earth_sun.getAcceleration(n-1);
        R = earth_sun.getDistance(n-1);
        printf("a = %.20f, r = %.3f \n", A[0], R[0]);
        earth_sun.writeToFile(filename + to_string(3));
        delete[] R; delete[] A;

        printf("v = %.3f\n", v4);
        planets[0].vy = v4;
        earth_sun.verletSolve(finalTime, n);
        A = earth_sun.getAcceleration(n-1);
        R = earth_sun.getDistance(n-1);
        printf("a = %.20f, r = %.3f \n", A[0], R[0]);
        earth_sun.writeToFile(filename + to_string(4));
        delete[] R; delete[] A;
    }
   earth_sun.destroy();
   delete[] planets;
}
