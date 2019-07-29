#ifndef MATRICES_H_
#define MATRICES_H_

#include <math.h>
#include <stdlib.h>

typedef struct SYMMAT33 {
    float data[6];
} SYMMAT33;

typedef struct MAT33 {
    float data[9];
} MAT33;

/****************************************************************************
 * Matrix / Vector routines
 ****************************************************************************/

SYMMAT33 *SYMMAT33_make(float a00, float a01, float a02,
			float a10, float a11, float a12,
			float a20, float a21, float a22);

MAT33 *MAT33_make(float a00, float a01, float a02,
		  float a10, float a11, float a12,
		  float a20, float a21, float a22);

void MAT33_assign(float a00, float a01, float a02,
		  float a10, float a11, float a12,
		  float a20, float a21, float a22, MAT33 *dest);

SYMMAT33 *MAT33_to_SYMMAT33(MAT33 *src, SYMMAT33 *dest);
SYMMAT33 *SYMMAT33_mult(SYMMAT33 *A, SYMMAT33 *B, SYMMAT33 *C);
MAT33 *MAT33_mult(MAT33 *A, MAT33 *B, MAT33 *C); 
float *VEC_crossp(float *a, float *b, float *c);
float *VEC_make_unit(float *v);


#endif /* MATRICES_H */
