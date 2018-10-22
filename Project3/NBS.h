#include <unordered_map>
#pragma once
#include <string>
#include "utils.h"

using namespace std;

struct MassObject {
	std::string name;
	double mass;
	// position:
	double x;
	double y;
	double z;
	//velocity:
	double vx;
	double vy;
	double vz;
};

class NBS {
protected:
	int m_n = 2;
	int m_m = 1;
	double pi = 3.14159265359;
	double **pos_x, **pos_y, **pos_z, **vel_x, **vel_y, **vel_z, **r;
	double *rel_pos_x, *rel_pos_y;
	double x0, y0, xn, yn;
	MassObject *massObjects;

	void distance(int);
	virtual void acceleration(int, int, double*, double*, double*) = 0;
	void createNBS();
	void createNBSrel();
	void setInit();
	void deleteNBS();

public:
	NBS(MassObject*, int); // constructor
	~NBS();
	void eulerSolve(double, int);
	void verletSolve(double, int);
	void verletSolveRel(double, int, int);
	void writeToFile(std::string);
	double timeEulerSolve(double, int);
	double timeVerletSolve(double, int);
	double* getAcceleration(int);
	double* getDistance(int);
};

class SolarSystem : public NBS {
protected:
	double pi = 3.14159265359;
	double fourPiPi = 4*pi*pi;
    double SM = 2e30;
	double m_centerMass = 0;
	double m_beta = 3;
	void acceleration(int, int, double*, double*, double*);
public:
	using NBS::NBS;
	void setCenterMass(double);
	void setBeta(double);
	void perihelionPrecession(int);
};

class SolarSystemRelativistic : public SolarSystem {
	void acceleration(int, int, double*, double*, double*);
	double c = 63239.7263; // Speed of light AU/yr
	double cc = c*c;
	double l;
public:
	SolarSystemRelativistic(MassObject*, int, int);
};


