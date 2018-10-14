#include <unordered_map>

std::unordered_map<std::string, double> Earth {
    {"mass", 6*10e24},
    {"x0", 9.413801075750535e-01},
    {"y0", 3.379019986046322E-01},
    {"z0", -9.334104672733438E-05},
    {"vx0", -5.994522787486753E-03*360.0},
    {"vy0", 1.617377250092178E-02*360.0},
    {"vz0", -1.732657683299539E-07*360.0}
};

std::unordered_map<std::string, double> Jupiter {
    {"mass", 1.9*10e27},
    {"x0", -2.666952709077877},
    {"y0", -4.655671225645230},
    {"z0", 7.896515774211305E-02},
    {"vx0", 6.458958874387921E-03*360.0},
    {"vy0", -3.390642961368397E-03*360.0},
    {"vz0", -1.303431975919576E-04*360.0}
};

std::unordered_map<std::string, double> Mars {
    {"mass", 6.6*10e23},
    {"x0", 1.377524608498332},
    {"y0", -1.476563536087052E-01},
    {"z0", -3.712365888122099E-02},
    {"vx0", 2.090430895774286E-03*360.0},
    {"vy0", 1.510431174964668E-02*360.0},
    {"vz0", 2.651531936464447E-04*360.0}
};

std::unordered_map<std::string, double> Venus {
    {"mass", 4.9*10e24},
    {"x0", 7.066148066072473E-01},
    {"y0", 1.672271996665292E-01},
    {"z0", -3.866246273128392E-02},
    {"vx0", -4.546772256640751E-03*360.0},
    {"vy0", 1.963882295765339E-02*360.0},
    {"vz0", 5.315229360444450E-04*360.0}
};

std::unordered_map<std::string, double> Saturn {
    {"mass", 5.5*10e26},
    {"x0", 1.549159416633817},
    {"y0", -9.935197133379326},
    {"z0", 1.110808051816651E-01},
    {"vx0", 5.204592319579984E-03*360.0},
    {"vy0", 8.422049097583516E-04*360.0},
    {"vz0", -2.220199495831516E-04*360.0}
};

std::unordered_map<std::string, double> Mercury {
    {"mass", 3.3*10e23},
    {"x0", -1.534743586808411E-01},
    {"y0", -4.321686982270908E-01},
    {"z0", -2.191308023125118E-02},
    {"vx0", 2.090296530646647E-02*360.0},
    {"vy0", -7.871323434203551E-03*360.0},
    {"vz0", -2.561512430419179E-03*360.0}
};

std::unordered_map<std::string, double> Uranus {
    {"mass", 8.8*10e25},
    {"x0", 1.717591494590517E+01},
    {"y0", 9.997664143031971},
    {"z0", -1.853846122526302E-01},
    {"vx0", -2.007356242188686E-03*360.0},
    {"vy0", 3.215850240122884E-03*360.0},
    {"vz0", 3.800786256690271E-05*360.0}
};

std::unordered_map<std::string, double> Neptune {
    {"mass", 1.03*10e26},
    {"x0", 2.892029941220658E+01},
    {"y0", -7.722539840090450},
    {"z0", -5.074661608355547E-01},
    {"vx0", 7.890393808745537E-04*360.0},
    {"vy0", 3.051931545808817E-03*360.0},
    {"vz0", -8.068388455453793E-05*360.0}
};

std::unordered_map<std::string, double> Pluto {
    {"mass", 1.31*10e22},
    {"x0", 1.164390186496279E+01},
    {"y0", -3.157511878129099E+01},
    {"z0", 1.062859645894982E-02},
    {"vx0", 3.018604420452015E-03*360.0},
    {"vy0", 4.214379702145380E-04*360.0},
    {"vz0", -9.301537706126110E-04*360.0}
};
