#ifndef LSSQ_H_
#define LSSQ_H_

#include "matrices.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef struct LSSQ_ARENA {
    double **btranspose_b;
    double **btranspose_b_inv;
    SYMMAT33 **bmat_tmp;
    double *xmean;
    double *sigmax;
    int nrows;
    int ncols;
} LSSQ_ARENA;

LSSQ_ARENA *lssq_initialize(SYMMAT33 **b_matrices, int nmatrices, int nentries);
void lssq_free(LSSQ_ARENA *lssq);
double lssq_compute(LSSQ_ARENA *lssq, double *X, double *B);
int lssq_gauss_jordan_solve (double **A, double **I, double *x, double *b, int n);

float *lssq_vector(int nl, int nh);
void lssq_free_vector(float *v, int nl, int nh);
void lssq_jacobi(float **a, int n, float d[], float **v, int *nrot,
		 float *b, float *z);
void lssq_eigen_sort(float d[], float **v, int n);

#endif /* LSSQ_H_ */
