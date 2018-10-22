#include <string>
#include "NBS.h"
#include "planets_v2.h"

using namespace std;
using namespace StellarObjectsLibraryv2;

void earthSun();
void allPlanets();

int main() {

	//earthSun();
	allPlanets();
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
