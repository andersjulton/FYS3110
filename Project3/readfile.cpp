#include <iostream>
#include <string>
#include <fstream>

void readPlanetData(std::string *planets, double **values, std::string filename, int numPlanets);

double **createMatrix(int m, int n) {
	double **mat;
	mat = new double*[m];
	for(int i = 0; i < m; i++) {
		mat[i] = new double[n];
		for(int j = 0; j < n; j++) {
			mat[i][j] = 0.0;
		}
	}
	return mat;
}

int main() {
    int numPlanets = 9;
    std::string filename = "Initial_values.txt";
    double **values;
    values = createMatrix(numPlanets, 6);
    std::string *planets;
    planets = new std::string[numPlanets];
    readPlanetData(planets, values, filename, numPlanets);

    std::cout << values[0][5];

    return 0;
}

void readPlanetData(std::string *planets, double **values, std::string filename, int numPlanets) {
    std::ifstream in(filename);
    if(in.is_open()) {
        double value;
        std::string planet;
        getline(in, planet);
        for (int i = 0; i < numPlanets; i++) {
            in >> planet;
            planets[i] = planet;
        }
        getline(in, planet);
        for (int i = 0; i < numPlanets; i++) {
            for (int j = 0; j < 6; j++) {
                in >> value;
                values[i][j] = value;
            }
        }
    }
}
