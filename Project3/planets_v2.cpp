#include <string> // not sure if this is needed.
#include "massObject.h"

namespace StellarObjectsLibraryv2 {

    struct MassObjectv2 Sun = {"Sun", 2e30, -1.474798697566110E-04, 7.248922398410132E-03, -7.268118257496195E-05,
        -7.585746627696994E-06*360.0, 2.590655894257502E-06*360.0, 1.896113821987181E-07*360.0};

    struct MassObjectv2 Mercury = {"Mercury", 3.3e23, -1.534743586808411E-01, -4.321686982270908E-01, -2.191308023125118E-02,
         2.090296530646647E-02*360.0, -7.871323434203551E-03*360.0, -2.561512430419179E-03*360.0};

    struct MassObjectv2 Venus = {"Venus", 4.9e24, 7.066148066072473E-01, 1.672271996665292E-01, -3.866246273128392E-02,
        -4.546772256640751E-03*360.0, 1.963882295765339E-02*360.0, 5.315229360444450E-04*360.0};

    struct MassObjectv2 Earth = {"Earth", 6e24, 9.413801075750535e-01, 3.379019986046322E-01, -9.334104672733438E-05,
        -5.994522787486753E-03*360.0, 1.617377250092178E-02*360.0, -1.732657683299539E-07*360.0};

    struct MassObjectv2 Mars = {"Mars", 6.6e23, 1.377524608498332, -1.476563536087052E-01, -3.712365888122099E-02,
         2.090430895774286E-03*360.0, 1.510431174964668E-02*360.0, 2.651531936464447E-04*360.0};

    struct MassObjectv2 Jupiter = {"Jupiter", 1.9e27, -2.666952709077877, -4.655671225645230, 7.896515774211305E-02,
        6.458958874387921E-03*360.0, -3.390642961368397E-03*360.0, -1.303431975919576E-04*360.0};

    struct MassObjectv2 Saturn = {"Saturn", 5.5e26, 1.549159416633817, -9.935197133379326, 1.110808051816651E-01,
         5.204592319579984E-03*360.0, 8.422049097583516E-04*360.0, -2.220199495831516E-04*360.0};

    struct MassObjectv2 Uranus = {"Uranus", 8.8e25, 1.717591494590517E+01, 9.997664143031971, -1.853846122526302E-01,
         -2.007356242188686E-03*360.0, 3.215850240122884E-03*360.0, 3.800786256690271E-05*360.0};

    struct MassObjectv2 Neptune = {"Neptune", 1.03e26, 2.892029941220658E+01, -7.722539840090450, -5.074661608355547E-01,
         7.890393808745537E-04*360.0, 3.051931545808817E-03*360.0, -8.068388455453793E-05*360.0};

    struct MassObjectv2 Pluto = {"Pluto", 1.31e22, 1.164390186496279E+01, -3.157511878129099E+01, 1.062859645894982E-02,
         3.018604420452015E-03*360.0, 4.214379702145380E-04*360.0, -9.301537706126110E-04*360.0};
}
