#pragma once

double *bisect(double *off, double *dia, int n);

int isolate(double xmin, double xmax, double *off, double *dia, double *eig, int count, int n);

double extract(double xmin, double xmax, double *off, double *dia, int n);

double *getInterval(double *off, double *dia, int n);

int getSignChange(double *a, double *b, int n, double lambda);
