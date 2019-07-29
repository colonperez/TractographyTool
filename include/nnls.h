#ifndef __NNLS_H
#define __NNLS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct NNLS_ARENA {
    double *A;
    double *A_orig;
    int n;
    int m;
    double *x;
    double *w;
    int mode;
    double rnorm;
    int *index;
    double *zz;
} NNLS_ARENA;

/****************************************************************************
 * NNLS routines
 ****************************************************************************/

NNLS_ARENA *nnls_initialize(double *A, int m, int n);
void nnls_compute(NNLS_ARENA *nnls, double *A, double *b, int copy_a);

int nnls_c(double* a, const int* mda, const int* m, const int* n, double* b, 
	   double* x, double* rnorm, double* w, double* zz, int* index, 
	   int* mode);


#endif /* __NNLS_H */