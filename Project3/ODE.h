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

class ODE {
protected:
	int m_n = 2;
	int m_m = 1;
	double pi = 3.14159265359;
	double **pos_x, **pos_y, **pos_z, **vel_x, **vel_y, **vel_z, **r;
	MassObject *massObjects;

	void distance(int);
	virtual void acceleration(int, int, double*, double*, double*) = 0;
	void createODE();
	void setInit();
	void deleteODE();
public:
	ODE(MassObject*, int); // constructor
	~ODE();
	void eulerSolve(double, int);
	void integrateVerlet(double, int);
	void writeToFile(std::string);
	double timeEulerSolve(double, int);
	double timeVerletSolve(double, int);
};

class SolarSystem : public ODE {
protected:
	void acceleration(int, int, double*, double*, double*);
public:
	using ODE::ODE;
};
