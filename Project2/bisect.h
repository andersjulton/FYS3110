#pragma once

double *Gershgorin(double *off, double *dia, int n);

int isolate(double xmin, double xmax, double *off, double *dia, double *eig, int count, int n);

double bisect(double xmin, double xmax, double *off, double *dia, int n);

double *getInterval(double *off, double *dia, int n);

int getSignChange(double *a, double *b, int n, double lambda);

double *getInterval(double *a, double *b, int n);
