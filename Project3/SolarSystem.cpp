#include "ODE.h"
#include "utils.h"


void SolarSystem::acceleration(int i, int j, double *ax, double *ay, double *az) {
    double fourPiPi = 4*pi*pi;
    double SM = 2e30;  // what is this?
    double  denum;;
	*ax = 0;
	*ay = 0;
	*az = 0;
	for (int k = 0; k < m_m; k++) {
		if (r[j][k] == 0) {
			continue;   // why? This means it only considers the earlier planets...
		}
		denum = SM*pow(r[j][k], 3);
		*ax -= fourPiPi*(pos_x[j][i] - pos_x[k][i])*massObjects[k].mass/denum;
		*ay -= fourPiPi*(pos_y[j][i] - pos_y[k][i])*massObjects[k].mass/denum;
		*az -= fourPiPi*(pos_z[j][i] - pos_z[k][i])*massObjects[k].mass/denum;
	}
}