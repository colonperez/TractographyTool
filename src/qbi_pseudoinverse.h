#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

gsl_matrix* moore_penrose_pinv(gsl_matrix *A, double rcond);