
int jacobi(double **A, double **R, int n);

int* getMaxInRow(double **A, int n);

double maxOffDiag(double **A, int *indexOfMax, int *k, int *l, int n);

void updateMaxInRow(double **A, int *indexOfMax, int k, int l, int n);

void rotate(double **A, double **R, int k, int l, int n);

