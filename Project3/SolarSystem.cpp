#include "NBS.h"
#include "utils.h"
#include <cmath>
#include <iostream>


void SolarSystem::acceleration(int i, int j, double *ax, double *ay, double *az) {
    double denum;
    double x = pos_x[j][i];
    double y = pos_y[j][i];
    double z = pos_z[j][i];
    double this_r = sqrt(x*x + y*y + z*z);
    denum = G*m_centerMass/pow(this_r, m_beta);
	*ax = -x*denum;
	*ay = -y*denum;
	*az = -z*denum;
	for (int k = 0; k < m_m; k++) {
		if (j == k) {
			continue;
		}
		denum = G*massObjects[k].mass/(pow(r[j][k], m_beta));
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


SolarSystemRelativistic::SolarSystemRelativistic(MassObject* initValue, int m) 
	: SolarSystem(initValue, m) {
	l = createVector(m, 0);
	double lx, ly, lz;
	MassObject planet;
	for (int i = 0; i < m; i++) {
		planet = massObjects[i];
		lx = planet.y*planet.vz - planet.z*planet.vy;
		ly = planet.x*planet.vz - planet.z*planet.vx;
		lz = planet.x*planet.vy - planet.y*planet.vx;
		l[i] = sqrt(lx*lx + ly*ly + lz*lz);
	}
}

void SolarSystemRelativistic::acceleration(int i, int j, double *ax, double *ay, double *az) {
    double  denum, rel;
    double threell_cc = 3*l[j]*l[j]/cc;
    double x = pos_x[j][i];
    double y = pos_y[j][i];
    double z = pos_z[j][i];
    double this_r = sqrt(x*x + y*y + z*z);
    rel = threell_cc/(this_r*this_r);
    denum = G*m_centerMass/pow(this_r, m_beta);
	*ax = -x*denum*(1 + rel);
	*ay = -y*denum*(1 + rel);
	*az = -z*denum*(1 + rel);
	for (int k = 0; k < m_m; k++) {
		if (j == k) {
			continue;
		}
		this_r = r[j][k];
		denum = G*massObjects[k].mass/(pow(this_r, m_beta));
		rel = threell_cc/(this_r*this_r);
		*ax -= (x - pos_x[k][i])*denum*(1 + rel);
		*ay -= (y - pos_y[k][i])*denum*(1 + rel);
		*az -= (z - pos_z[k][i])*denum*(1 + rel);
	}
}

void SolarSystemRelativistic::destroy() {
	SolarSystem::destroy();
	delete[] l;
}
