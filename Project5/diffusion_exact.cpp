#define _USE_MATH_DEFINES

#include "utils.h"
#include <cmath>

double u(double x, double t, double L, int n) {
    double val = 0;
    for(int i = 1; i < n; i++) {
        val += pow(-1, i+1)*2*L/(M_PI*n)*sin(i*M_PI*x/L)*exp(-pow(i*M_PI, 2)*t)
    }
    return x - val;
}

double *exact(double *v, double x0, double x1, double m, double t, double L, int n) {
    double dx = (x1 - x0)/m;
    double x;
    for(int i = 0; i < (m+1); i++) {
        x = x0 + i*dx
        v[i] = u(x, t, L, n);
    }
}


// NOT DONE
