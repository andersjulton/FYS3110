#include "NBS.h"
#include "utils.h"
#include <cmath>
#include <iostream>


void SolarSystem::acceleration(int i, int j, double *ax, double *ay, double *az) {
    double denum;
    double con = fourPiPi/SM;
    double x = pos_x[j][i];
    double y = pos_y[j][i];
    double z = pos_z[j][i];
    double this_r = sqrt(x*x + y*y + z*z);
    denum = fourPiPi*m_centerMass/(SM*pow(this_r, m_beta));
	*ax = -x*denum;
	*ay = -y*denum;
	*az = -z*denum;
	for (int k = 0; k < m_m; k++) {
		if (j == k) {
			continue;
		}
		denum = con*massObjects[k].mass/(pow(r[j][k], m_beta));
		*ax -= (x - pos_x[k][i])*denum;
		*ay -= (y - pos_y[k][i])*denum;
		*az -= (z - pos_z[k][i])*denum;
	}
}

void SolarSystem::setCenterMass(double centerMass) {
	m_centerMass = centerMass;
}

void SolarSystem::setBeta(double beta) {
	m_beta = beta + 1;
}

void SolarSystem::perihelionPrecession(int years) {
	
	double r;
	double rmin;
	int index = 0;

	rmin = sqrt(pow(rel_pos_x[0], 2) + pow(rel_pos_y[0], 2));
	for (int i = 0; i < m_n/years; i++) {
		r = sqrt(pow(rel_pos_x[i], 2) + pow(rel_pos_y[i], 2));
		if (r < rmin) {
			rmin = r;
			index = i;
		}
	}
	std::cout << "Position of perihelion is: (" << rel_pos_x[index] << "," << rel_pos_y[index] << ")\n";
	std::cout << "Distance from sun: " << sqrt(pow(rel_pos_x[index], 2) + pow(rel_pos_y[index], 2)) << "\n";
	std::cout << "Argument of perihelion is: " << atan2(rel_pos_y[index], rel_pos_x[index])*206264.806 << " arcseconds/century.\n";

}


SolarSystemRelativistic::SolarSystemRelativistic(MassObject* initValue, int m, int index) 
	: SolarSystem(initValue, m) {
	double lx, ly, lz;
	MassObject planet;
	planet = massObjects[index];
	lx = planet.y*planet.vz - planet.z*planet.vy;
	ly = planet.x*planet.vz - planet.z*planet.vx;
	lz = planet.x*planet.vy - planet.y*planet.vx;
	l = sqrt(lx*lx + ly*ly + lz*lz);
	
}

void SolarSystemRelativistic::acceleration(int i, int j, double *ax, double *ay, double *az) {
    double denum, rel;
    double const1 = fourPiPi/SM;

    double const2 = 3*l*l/cc;
 
    double this_r = sqrt(x0*x0 + y0*y0);
    rel = const2/(this_r*this_r);
    denum = fourPiPi*m_centerMass/(SM*pow(this_r, m_beta));
	*ax = -x0*denum*(1 + rel);
	*ay = -y0*denum*(1 + rel);

}
