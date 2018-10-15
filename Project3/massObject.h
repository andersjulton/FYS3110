#include <unordered_map>
#include <string>

using namespace std;

class MassObject {
	public:
		double mass;
		// position:
		double x;
		double y;
		double z;
		//velocity:
		double vx;
		double vy;
		double vz;
		// constructor
    	MassObject(unordered_map<std::string, double>);
};