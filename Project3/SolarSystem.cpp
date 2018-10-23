#include "NBS.h"
#include "utils.h"
#include <cmath>
#include <iostream>


void SolarSystem::acceleration(int i, int j, double *ax, double *ay, double *az) {
    double denum;
    double x = pos_x[j][i];
    double y = pos_y[j][i];
    double z = pos_z[j][i];
    denum = G*m_centerMass/pow(r[j][j], m_beta);
	*ax = -x*denum;
	*ay = -y*denum;
	*az = -z*denum;
	for (int k = 0; k < m_m; k++) {
		if (j == k) {
			continue;
		}
		denum = G*massObjects[k].mass/(pow(r[j][k], m_beta));
		*ax += (pos_x[k][i] - x)*denum;
		*ay += (pos_y[k][i] - y)*denum;
		*az += (pos_z[k][i] - z)*denum;
	}
}

double* SolarSystem::kineticEnergy() {
	int h = (int) m_n/(m_finalTime + 1);
	int i;
	double *K = createVector(0, ((int) m_finalTime + 1));
	for (int index = 0; index < ((int) m_finalTime + 1); index++) {
		i = h*index;
		distance(i);
		for (int j = 0; j < m_m; j++) {
			K[index] += (vel_x[j][i]*vel_x[j][i] + vel_y[j][i]*vel_y[j][i] + vel_z[j][i]*vel_z[j][i])*massObjects[j].mass;
		}
		K[index] = K[index]*0.5;
	} 
	return K;
}

/*
// return total K, U and E for every year
double* SolarSystem::conservationEnergy() {  
	double *U = createVector(0, n);
	double *E = createVector(0, n);
	distance(0);

	double E_0 = K*0.5 - U*G;
	double E_min = E_0;
	double E_max = E_0;
	double E;
	printf("E0 = %f\n", E_0);
	for (int i = 1; i < m_n; i++) {
		K = 0;
		U = 0;
		distance(i);
		for (int j = 0; j < m_m; j++) {
			K += (vel_x[j][i]*vel_x[j][i] + vel_y[j][i]*vel_y[j][i] + vel_z[j][i]*vel_z[j][i])*massObjects[j].mass;
			for (int k = 0; k < m_m; k++) {
				if (j == k) {
					continue;
				}
				U += massObjects[j].mass*massObjects[k].mass/r[j][k];
			}
		} 
		E = K*0.5 - U*G;
	}
	printf("En = %f\n", E);
	printf("E_min = %f, E_max = %f\n", E_min, E_max);
	return (E_max - E_min)/E_0;
}
*/


double SolarSystem::conservationEnergy2() {
	double K = 0;  
	double U = 0;
	distance(0);
	for (int i = 1; i < m_m; i++) {
		K += (vel_x[i][0]*vel_x[i][0] + vel_y[i][0]*vel_y[i][0] + vel_z[i][0]*vel_z[i][0])*massObjects[i].mass;
		U += massObjects[i].mass/r[0][i];
	} 
	double E_0 = K*0.5 - U*fourPiPi;
	double E_min = E_0;
	double E_max = E_0;
	printf("E0 = %f\n", E_0); 
	double E;
	for (int i = 1; i < m_n; i++) {
		K = 0;
		U = 0;
		distance(i);
		for (int j = 1; j < m_m; j++) {
			K += pow((vel_x[j][i] - vel_x[0][i]), 2)*massObjects[j].mass;
			K += pow((vel_y[j][i] - vel_y[0][i]), 2)*massObjects[j].mass;
			K += pow((vel_z[j][i] - vel_z[0][i]), 2)*massObjects[j].mass;
			U += massObjects[j].mass/r[0][j];
			for (int k = 1; k < m_m; k++) {
				if (j == k) {
					continue;
				}
				U += massObjects[j].mass*massObjects[k].mass/r[j][k];
			}
		}
		E = K*0.5 - U*fourPiPi;
		//printf("%f\n", E);
		if (E < E_min) {
			E_min = E;
		} else if (E > E_max) {
			E_max = E;
		}
	}
	printf("En = %f\n", E);
	printf("E_min = %f, E_max = %f\n", E_min, E_max);
	return (E_max - E_min);
}

// from center of mass
double SolarSystem::conservationEnergy3() {
	double totalMass = m_centerMass;
	for (int i = 0; i < m_m; i++) {
		totalMass += massObjects[i].mass;
	}
	double K = 0;  
	double U = 0;
	distance(0);
	for (int i = 0; i < m_m; i++) {
		K += (vel_x[i][0]*vel_x[i][0] + vel_y[i][0]*vel_y[i][0] + vel_z[i][0]*vel_z[i][0])*massObjects[i].mass;
		U += totalMass*massObjects[i].mass/r[i][i];
	} 
	double E_0 = K*0.5 - U*G;
	double E_min = E_0;
	double E_max = E_0;
	printf("E0 = %f\n", E_0); 
	double E;
	for (int i = 1; i < m_n; i++) {
		K = 0;
		U = 0;
		distance(i);
		for (int j = 0; j < m_m; j++) {
			K += (vel_x[j][i]*vel_x[j][i] + vel_y[j][i]*vel_y[j][i] + vel_z[j][i]*vel_z[j][i])*massObjects[j].mass;
			U += totalMass*massObjects[j].mass/r[j][j];
		}
		E = K*0.5 - U*G;
		//printf("%f\n", E);
		if (E < E_min) {
			E_min = E;
		} else if (E > E_max) {
			E_max = E;
		}
	}
	printf("En = %f\n", E);
	printf("E_min = %f, E_max = %f\n", E_min/E_0, E_max/E_0);
	return (E_max - E_min)/E_0;
}

void SolarSystem::setCenterMass(double centerMass) {
	m_centerMass = centerMass;
}

void SolarSystem::setBeta(double beta) {
	m_beta = beta + 1;
}

void SolarSystemRelativistic::perihelionPrecession(int index, double finalTime, int n, int years) {
	double *per;
	per = verletSolveRel2D(index, finalTime, n, years);
	double r = sqrt(pow(per[0], 2) + pow(per[1], 2));

	std::cout << "Position of perihelion is: (" << per[0] << "," << per[1] << ")\n";
	std::cout << "Distance from sun: " << r << "\n";
	std::cout << "Argument of perihelion is: " << atan2(per[1], per[0])*206264.806 << " arcseconds/century.\n";
	delete[] per;
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

double *SolarSystemRelativistic::verletSolveRel2D(int index, double finalTime, int n, int years) {
	double x0, xn, y0, yn, vx0, vy0, vxn, vyn, a, an;
	double r, rn, rmin;
	double *perihelion;
	double ll = l[index]*l[index];
	perihelion = createVector(0, 2);

	x0 = massObjects[index].x;
	y0 = massObjects[index].y;
	vx0 = massObjects[index].vx;
	vy0 = massObjects[index].vy;
	distance(0);

	double h = finalTime / (n + 1);
	double hh = h * h;
	int count = 0;
	r = sqrt(pow(x0, 2) + pow(y0, 2));
	a = -fourPiPi * (1.0 + 3.0*ll / (r*r*cc)) / pow(r, 3);

	for (int i = 0; i < n - n / 100; i++) {
		xn = x0 + h*vx0 + (hh / 2.0)*a*x0;
		yn = y0 + h*vy0 + (hh / 2.0)*a*y0;

		rn = sqrt(xn*xn + yn*yn);
		an = -fourPiPi*(1.0 + 3.0*ll/(rn*rn*cc))/pow(rn, 3);

		vxn = vx0 + (h / 2.0)*(an*xn + a*x0);
		vyn = vy0 + (h / 2.0)*(an*yn + a*y0);

		r = rn;
		a = an;
		x0 = xn;
		y0 = yn;
		vx0 = vxn;
		vy0 = vyn;
	}

	rmin = sqrt(pow(x0, 2) + pow(y0, 2));
	perihelion[0] = x0;
	perihelion[1] = y0;

	for (int i = 0; i < n / 100; i++) {
		xn = x0 + h * vx0 + (hh / 2.0)*a*x0;
		yn = y0 + h * vy0 + (hh / 2.0)*a*y0;

		r = sqrt(pow(xn, 2) + pow(yn, 2));
		if (r < rmin) {
			rmin = r;
			perihelion[0] = xn;
			perihelion[1] = yn;
		}

		rn = sqrt(xn*xn + yn*yn);
		an = -fourPiPi * (1.0 + 3.0*ll/(rn*rn*cc))/pow(rn, 3);

		vxn = vx0 + (h / 2.0)*(an*xn + a * x0);
		vyn = vy0 + (h / 2.0)*(an*yn + a * y0);

		r = rn;
		a = an;
		x0 = xn;
		y0 = yn;
		vx0 = vxn;
		vy0 = vyn;
	}
	return perihelion;
}
