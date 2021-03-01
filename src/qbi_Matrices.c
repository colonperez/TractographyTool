#include"qbi_Matrices.h"
#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include"qbi_pseudoinverse.h"

double mulPoints(double* rowM, double* colM, int length){
	//printf("mulPoints\n");
	double out = 0;
	for(int i = 0; i < length; ++i){
		out += rowM[i] * colM[i];
		//printf("\t\tmulPoints i = %d\n", i);
	}
	/*if(isnan(out)){
		printf("Expression is nan. Inputs:\t");
		for(int i = 0; i < length; ++i){
			printf("%f\t",rowM[i]);
		}
		printf(" and ");
		for(int i = 0; i < length; ++i){
			printf("%f\t",colM[i]);
		}
		printf("\n");
	}*/
	return out;
}

void matMulPoints(double p1[3], double p2[3], double out[3][3]){

	out[0][0] = p1[0] * p2[0];
	out[0][1] = p1[0] * p2[1];
	out[0][2] = p1[0] * p2[2];

	out[1][0] = p1[1] * p2[0];
	out[1][1] = p1[1] * p2[1];
	out[1][2] = p1[1] * p2[2];

	out[2][0] = p1[2] * p2[0];
	out[2][1] = p1[2] * p2[1];
	out[2][2] = p1[2] * p2[2];
}

void pseudoInverse(double**H, int rowDim, int colDim, double**out){
	//Converts two dimensional array to libgsl matrix, then calls moore_penrose_pinv and converts the output back to a two dimensional array

	gsl_matrix* A = gsl_matrix_alloc(rowDim,colDim);

	for(int row = 0; row < rowDim; ++row){
		for(int col = 0; col < colDim; ++col){
			gsl_matrix_set(A,row,col,H[row][col]);
		}
	}

	gsl_matrix* o = moore_penrose_pinv(A,1E-15);

	gsl_matrix_free(A);

	for(int row = 0; row < colDim; ++row){
		for(int col = 0; col < rowDim; ++col){
			out[row][col] = gsl_matrix_get(o,row,col);
		}
	}

	gsl_matrix_free(o);

}

double mean(double* ar, int n){
	double m = 0;
	for(int i = 0; i < n; ++i){
		m += ar[i];
	}
	m /= n;
	return m;
}

double rms(double* ar, int n){
	double sum = 0;

	for (int i = 0; i < n; i++)
		sum += pow(ar[i], 2);

	return sqrt(sum / n);
}

double std(double* ar, int n){
	double m = mean(ar,n);
	double sum;
	for(int i = 0; i < n; ++i){
		sum+= ((ar[i] - m)*(ar[i] - m));
	}
	sum /= (n-1);
	return sqrt(sum);
}