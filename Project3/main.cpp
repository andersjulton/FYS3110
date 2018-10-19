#include <unordered_map>
#include <string>
#include "ODE.h"
#include "planets_v2.h"
using namespace std;


int main() {
	int m = 10;
    MassObject *planets = new MassObject[m];

    int n = 100000;
    double finalTime = 250.0;

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
    ths.writeToFile();

    delete[] planets;
	return 0;
}
