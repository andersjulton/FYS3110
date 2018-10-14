#include <unordered_map>
#include <iostream>
#include <string>

class Planet(std::unordered_map<std::string, double> data) {
private:
    double mass = data["mass"];
    double x0 = data["x0"];
    double y0 = data["y0"];
    double z0 = data["z0"];
    double vx0 = data["vx0"];
    double vy0 = data["vy0"];
    double vz0 = data["vz0"];
};
