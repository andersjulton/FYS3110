#include <string>
#include "NBS.h"
#include "planets_v2.h"

using namespace std;
using namespace StellarObjectsLibraryv2;

void earthSun();
void allPlanets();
void escapeVelocity();

int main() {

	escapeVelocity();
	return 0;
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
    double v_e = 8.88576587631731840133;
    double v = 8.875;
    MassObject planet = {"Earth", 6e24, 1, 0, 0, 0, v, 0};
    planets[0] = planet;
    SolarSystem earth_sun(planets, m);
    earth_sun.setCenterMass(Sun.mass);
    finalTime = 3500.0;
    n = 1000000;
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

    printf("Analytical %.20f \n", sqrt(8)*pi);
    earth_sun.writeToFile("Earth_Sun_min");

    v = 8.89;
    planets[0].vy = v;
    finalTime = 3500.0;
    n = 1000000;
    SolarSystem earth_sun2(planets, m);
    earth_sun2.setCenterMass(Sun.mass);
    earth_sun2.verletSolve(finalTime, n);

    A = earth_sun2.getAcceleration(1);
    R = earth_sun2.getDistance(1);
    printf("a = %.20f, r = %.3f \n", A[0], R[0]);
    A = earth_sun2.getAcceleration(n/2);
    R = earth_sun2.getDistance(n/2);
    printf("a = %.20f, r = %.3f \n", A[0], R[0]);
    A = earth_sun2.getAcceleration(n-1);
    R = earth_sun2.getDistance(n-1);
    printf("a = %.20f, r = %.3f \n", A[0], R[0]);

    earth_sun2.writeToFile("Earth_Sun_max");

    delete[] planets;
    delete[] R;
    delete[] A;
}
