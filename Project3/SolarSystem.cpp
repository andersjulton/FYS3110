#include "ODE.h"
#include "utils.h"
#include <cmath>


void SolarSystem::acceleration(int i, int j, double *ax, double *ay, double *az) {
    double fourPiPi = 4*pi*pi;
    double SM = 2e30;
    double  denum;
    double x = pos_x[j][i];
    double y = pos_y[j][i];
    double z = pos_z[j][i];
    double this_r = sqrt(x*x + y*y + z*z);
    denum = fourPiPi*m_centerMass/SM*pow(this_r, 3);
	*ax = x*denum;
	*ay = y*denum;
	*az = z*denum;
	for (int k = 0; k < m_m; k++) {
		if (j == k) {
			continue;   
		}
		denum = fourPiPi*massObjects[k].mass/(SM*pow(r[j][k], 3));
		*ax -= (x - pos_x[k][i])*denum;
		*ay -= (y - pos_y[k][i])*denum;
		*az -= (z - pos_z[k][i])*denum;
	}
}

void SolarSystem::setCenterMass(double centerMass) {
	m_centerMass = centerMass;
}