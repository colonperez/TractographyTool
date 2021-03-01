
double mulPoints(double* rowM, double* colM, int length);

void matMulPoints(double p1[3], double p2[3], double out[3][3]);

void pseudoInverse(double ** mat, int rowDim, int colDim, double ** out);

double mean(double* ar, int n);
double rms(double* ar, int n);
double std(double* ar, int n);