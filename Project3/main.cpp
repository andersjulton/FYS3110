#include <unordered_map>
#include <string>
#include "NBS.h"
#include "planets_v2.h"

using namespace std;
using namespace StellarObjectsLibraryv2;

int main() {
	int m, n;
	double finalTime;
	string filename;
	MassObject *planets;

	// part a) Earth-Sun system, fixed mass center.
	m = 1;
    planets = new MassObject[m];
    // what should the initial velocity be here?
    MassObject Earth2D = {"2DEarth", 6e24, 1, 0, 0,
        -5.994522787486753E-03*360.0, 1.617377250092178E-02*360.0, 0};
    planets[0] = Earth2D;
    finalTime = 50.0;
    n = finalTime*1000;
    SolarSystem earth_sun(planets, m);
    earth_sun.setCenterMass(Sun.mass);
    filename = "earth_sun_euler";
    earth_sun.eulerSolve(finalTime, n);
    earth_sun.writeToFile(filename);
    filename = "earth_sun_verlet";
    earth_sun.verletSolve(finalTime, n);
    earth_sun.writeToFile(filename);
    delete[] planets;

    // all planets
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


	return 0;
}
