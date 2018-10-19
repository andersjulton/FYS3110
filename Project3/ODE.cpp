#include <unordered_map>
#include <string>
#include "ODE.h"
#include "utils.h"
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

ODE::ODE(MassObject* initValue, int new_m) {
	m = new_m;
	massObjects = initValue;
	createODE();
	setInit();
	r = createMatrix(m, m);
	distance(0);
}

void ODE::eulerSolve(int n) {
    double h = 1.0/(n+1);
    double r;
    double fourPiPi = 4*pi*pi;
    for (int i = 0; i < (n-1); i++) {
        r = sqrt(pos_x[1][i]*pos_x[1][i] + pos_y[1][i]*pos_y[1][i]);

        pos_x[1][i+1] = pos_x[1][i] + h*vel_x[1][i];
        pos_y[1][i+1] = pos_y[1][i] + h*vel_y[1][i];

        vel_x[1][i+1] = vel_x[1][i] - h*fourPiPi/(r*r*r)*pos_x[1][i];
        vel_y[1][i+1] = vel_y[1][i] - h*fourPiPi /(r*r*r)*pos_y[1][i];
    }
}

void ODE::integrateVerlet(double finalTime, int steps) {
	deleteODE();
	n = steps;
	createODE();
	setInit();
	distance(0);

	h = finalTime/(n+1);
    double hh = h*h;
    double ax, ay, az;
    double **A;
    A = createMatrix(m, 3);
    for (int i = 0; i < m; i++) {
    	acceleration(0, i, &ax, &ay, &az);
        A[i][0] = ax;
        A[i][1] = ay;
        A[i][2] = az;

    }
    for (int i = 0; i < n-1; i++) {
        for (int j = 0; j < m; j++) {
            pos_x[j][i+1] = pos_x[j][i] + h*vel_x[j][i] + (hh/2.0)*A[j][0];
            pos_y[j][i+1] = pos_y[j][i] + h*vel_y[j][i] + (hh/2.0)*A[j][1];
            pos_z[j][i+1] = pos_z[j][i] + h*vel_z[j][i] + (hh/2.0)*A[j][2];
            acceleration((i+1), j, &ax, &ay, &az);
            vel_x[j][i+1] = vel_x[j][i] + (h/2.0)*(A[j][0] + ax);
            vel_y[j][i+1] = vel_y[j][i] + (h/2.0)*(A[j][1] + ay);
            vel_z[j][i+1] = vel_z[j][i] + (h/2.0)*(A[j][2] + az);

            A[j][0] = ax;
            A[j][1] = ay;
            A[j][2] = az;
        }
        distance(i+1);
    }
    deleteMatrix(A, m);
}

void ODE::deleteODE() {
	deleteMatrix(pos_x, m);
	deleteMatrix(pos_y, m);
	deleteMatrix(pos_z, m);
	deleteMatrix(vel_x, m);
	deleteMatrix(vel_y, m);
	deleteMatrix(vel_z, m);
}

void ODE::writeToFile() {
	doubleMatrixToFile(pos_x, n, m, "pos_xTEST");
    doubleMatrixToFile(pos_y, n, m, "pos_yTEST");
    doubleMatrixToFile(pos_z, n, m, "pos_zTEST");
}

void ODE::createODE() {
	pos_x = createMatrix(m, n);
    pos_y = createMatrix(m, n);
    pos_z = createMatrix(m, n);
    vel_x = createMatrix(m, n);
    vel_y = createMatrix(m, n);
    vel_z = createMatrix(m, n);
}

void ODE::setInit(){
	for (int i = 0; i < m; i++) {
		pos_x[i][0] = massObjects[i].x;
		pos_y[i][0] = massObjects[i].y;
		pos_z[i][0] = massObjects[i].z;
		vel_x[i][0] = massObjects[i].vx;
		vel_y[i][0] = massObjects[i].vy;
		vel_z[i][0] = massObjects[i].vz;
	}
}

void ODE::distance(int iter) {
	double temp;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < i; j++) {
			temp = pow(pos_x[i][iter] - pos_x[j][iter], 2.0);
			temp += pow(pos_y[i][iter] - pos_y[j][iter], 2.0);
			temp += pow(pos_z[i][iter] - pos_z[j][iter], 2.0);
			temp = sqrt(temp);
			r[i][j] = temp;
			r[j][i] = temp;		
		}
	}
}

ODE::~ODE() {
	deleteODE();
	deleteMatrix(r, m);
}

void SolarSystem::acceleration(int i, int j, double *ax, double *ay, double *az) {
    double fourPiPi = 4*pi*pi;
    double SM = 2e30;  // what is this?
    double  denum;;
	*ax = 0;
	*ay = 0;
	*az = 0;
	for (int k = 0; k < m; k++) {
		if (r[j][k] == 0) {
			continue;   // why? This means it only considers the earlier planets...
		}
		denum = SM*pow(r[j][k], 3);
		*ax -= fourPiPi*(pos_x[j][i] - pos_x[k][i])*massObjects[k].mass/denum;
		*ay -= fourPiPi*(pos_y[j][i] - pos_y[k][i])*massObjects[k].mass/denum;
		*az -= fourPiPi*(pos_z[j][i] - pos_z[k][i])*massObjects[k].mass/denum;
	}
}