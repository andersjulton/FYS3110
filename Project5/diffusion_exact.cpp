#include "utils.h"
#include <cmath>
#define _USE_MATH_DEFINES

double u(double x, double t, int n) {
    double val = 0;
    for(int i = 1; i < n; i++) {
        val += pow(-1, i)*2/(M_PI*i)*sin(i*M_PI*x)*exp(-pow(i*M_PI, 2)*t);
    }
    return (x + val);
}

void analytic1D(double *v, double t, double dx, int n) {
    for(int i = 0; i < n; i++) {
        v[i] = u((i*dx), t, n);
    }
}


// NOT DONE
