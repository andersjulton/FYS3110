#include <unordered_map>
#include <string>
#include "massObject.h"

using namespace std;

MassObject::MassObject(unordered_map<std::string, double> map) {
	mass = map["mass"];
	x = map["x0"];
	y = map["y0"];
	z = map["y0"];
	vx = map["vx0"];
	vy = map["vy0"];
	vz = map["vz0"];
}