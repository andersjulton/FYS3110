#define _USE_MATH_DEFINES

#include "utils.h"
#include <cmath>

double u(double x, double t, int n) {
    double val = 0;
    for(int i = 1; i < n; i++) {
        val += pow(-1, i)*2/(M_PI*i)*sin(i*M_PI*x)*exp(-pow(i*M_PI, 2)*t)
    }
    return x + val;
}

void *analytic1D(double *v, double t, int n) {
    double dx = 1/n;
    double x;
    for(int i = 0; i < n; i++) {
        x = i*dx
        v[i] = u(x, t, n);
    }
}


// NOT DONE
