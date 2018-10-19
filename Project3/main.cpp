#include <unordered_map>
#include <string>
#include "ODE.h"
#include "planets_v2.h"
using namespace std;


int main() {
	int m = 2;
    MassObject *planets = new MassObject[m];

    int n = 100000;
    double finalTime = 250.0;

    string filename = "euler";

    planets[0] = StellarObjectsLibraryv2::Sun;
    planets[1] = StellarObjectsLibraryv2::Earth;

    SolarSystem hrrr(planets, m);
    hrrr.eulerSolve(finalTime, n);
    hrrr.writeToFile(filename);

    delete[] planets;

    m = 10;
    planets = new MassObject[m];
    filename = "pos";


    planets[0] = StellarObjectsLibraryv2::Sun;
    planets[1] = StellarObjectsLibraryv2::Mercury;
    planets[2] = StellarObjectsLibraryv2::Venus;
    planets[3] = StellarObjectsLibraryv2::Earth;
    planets[4] = StellarObjectsLibraryv2::Mars;
    planets[5] = StellarObjectsLibraryv2::Jupiter;
    planets[6] = StellarObjectsLibraryv2::Saturn;
    planets[7] = StellarObjectsLibraryv2::Uranus;
    planets[8] = StellarObjectsLibraryv2::Neptune;
    planets[9] = StellarObjectsLibraryv2::Pluto;

    SolarSystem ths(planets, m);
    ths.integrateVerlet(finalTime, n);
    ths.writeToFile(filename);

    double time = ths.timeVerletSolve(finalTime, n);
    printf("verlet = %f\n", time);
    time = ths.timeEulerSolve(finalTime, n);
    printf("euler = %f\n", time);

    delete[] planets;
	return 0;
}