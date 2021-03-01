#include"qbi_ConvolutionAndPrediction.h"
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"qbi_Matrices.h"
#include"qbi_pseudoinverse.h"

#define pi 3.14159265358979323846

void getDiffusionSignal(double **q, double **v, int qLen, int vLen,double**out){
	//INPUT: Two arrays of points, and their respective lengths
	//Out: Two dimensional m x p, equal to phi(cos-1(abs(Qt*V)))

	//NOTE: Also used in finding diffusion signal estimate

	double* tempQ = malloc(sizeof(double)*3);
	double* tempV = malloc(sizeof(double)*3);

	for(int i = 0; i < qLen; ++i){

		tempQ[0] = q[0][i];
		tempQ[1] = q[1][i];
		tempQ[2] = q[2][i];

		for(int j = 0; j < vLen; ++j){

			tempV[0] = v[0][j]; 
			tempV[1] = v[1][j]; 
			tempV[2] = v[2][j];

			double z = getKernel(tempQ,tempV);

			out[i][j] = z;
		}
	}

	free(tempQ);
	free(tempV);

}

void getReconstructionMatrix(double ** G, double ** H, int k, int n, int p, int m, double ** out) {
	//Result: (In X 1kT)GH+
		//X is the matrix direct product
	//the summation should be implemented by repartitioning the (kn) x m matrix GH+ into a k x n x m array and then summing over the first dimension.

	//Dimensions of G: (kn) x p
	//Dimensions of H: m x p
	//Dimensions of H+: p x m
	//Dimensions of Summation Matrix: n x kn
	//Dimensions of GH+: kn x m

	//Dimensions of out:
		//n x m


	double ** mul = malloc(sizeof(double*) * k * n);
	for(int i = 0; i < k * n; ++ i){
		mul[i] = malloc(sizeof(double) * m);
	}

	double ** HPlus = malloc(sizeof(double*) * p);
	for(int i = 0; i < p; ++i){
			HPlus[i] = malloc(sizeof(double) * m);
		}

	double** summation = malloc(sizeof(double*)*n);
	for(int i = 0; i < n; ++i){
		summation[i] = malloc(sizeof(double)*k*n);
	}

	int marker = 0;
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < k*n; ++j){
			if(j >= marker && j < (marker+k)){
				summation[i][j] = 1;
			}
			else{
				summation[i][j] = 0;
			}
		}
		marker += k;
	}

	pseudoInverse(H,m, p, HPlus);

	printf("H: \n");
	for(int i = 0; i < 5; ++i){
		for(int j = 0; j < 5; ++j){
			printf("%f,",H[i][j]);
		}
		printf("\n");
	}


	printf("H+: \n");
	for(int i = 0; i < 5; ++i){
        for(int j = 0; j < 5; ++j){
            printf("%f,",HPlus[i][j]);
        }
        printf("\n");
    }

    /*printf("G: \n");
    for(int i = 0; i < k*n; ++i){
    	for(int j = 0; j < p; ++j){
    		printf("%f ", G[i][j]);
    	}
    	printf("\n");
    }*/


	//Calculate GH+
	for(int row = 0; row < (k*n); row++){
		for(int col = 0; col < m; col++){
			double sum = 0;
			for(int mark = 0; mark < p; mark++){
				sum += G[row][mark] * HPlus[mark][col];
			}
			mul[row][col] = sum;
		}
	}

	printf("GH+: \n");
	for(int i = 0; i < 5;++i){
		for(int j = 0; j < 5; ++j){
			printf("%f,",mul[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	//Calculate (In X 1kT)GH+
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < m; ++j){
			double sum = 0;
			for(int _ = 0; _ < k*n; ++_){
				sum += summation[i][_]*mul[_][j];
			}
			out[i][j] = sum;
		}
	}

	for(int i = 0; i < (k * n); ++ i){
		free(mul[i]);
	}
	free(mul);

	for(int i = 0; i < p; ++i){
		free(HPlus[i]);
	}
	free(HPlus);

	for(int i = 0; i < n; ++i){
		free(summation[i]);
	}
	free(summation);
}


double calculateNormalizationConstant(double ** A, double * e, int n, int m){
	//Calculates the normalization constant Z = 1nTAe

	// 1nT = 1 x n vector of ones
	//A = Reconstruction matrix (use getReconstructionMatrix)
		//n x m
	//e = Diffusion Signal


	//Calculate 1nT*A
	//Resulting vector: 1 x m
	//(Sum each col)

	//printf("calculateNormalizationConstant: \n");

	double firstProduct[m];
	for(int i = 0; i < m; ++i){
		double sum = 0;
		for(int j = 0; j < n; ++j){
			sum += A[j][i];
		}
		firstProduct[i] = sum;
	}

	/*
	printf("First Matrix Product:\n");
	for(int i = 0; i < m; ++i) {
		printf("%lf\t",firstProduct[i]);
	}
	printf("\n");

	*/

	//Calculate firstProduct * e
	double out = 0;
	for(int i = 0; i < m; ++i) {
		out += e[i] * firstProduct[i];
	}

	return out;
}

void computeODF(double ** A, double * e, int n, int m, double * out) { //Normalization constant might not be necessary, since we are looking for peaks rather than values
	//Inputs: Reconstruction matrix A, diffusion signal e

	//Calculates the estimated ODF, (1/Z)Ae, with Z being a normalization constant

	//A: n x m
	//e: m x 1

	//Output: (1/Z)Ae, n x 1 vector

	//Normalization Constant
	double norm = 1 / calculateNormalizationConstant(A,e,n,m);

	//Calculate Ae and multiply each term by norm
	for(int i = 0; i < n; ++i){
		double count = 0;
		for(int j = 0; j < m; ++j){
			count += e[j] * A[i][j];
		}
		out[i] = count*norm;
		//printf("ODF[%d]: %f\t",i,count);
	}
	//printf("\n");

}

void getEquator(const int k, double ** C){

	//In:
	//	k sampling directions
	//	3 x k output matrix C

	//Operation:
	//	Takes k sampling directions and returns an equator of k equally spaced points

	//C must be 3 x k
	//To allocate space for c (use in main):

	for(int i = 0; i < k; ++i){
		C[0][i] = cos((2 * pi / k) * (i + 1));
		C[1][i] = sin((2 * pi / k) * (i + 1));
		C[2][i] = 0;
	}
}

void getRotationMat(double z[3],double u[3], double** out){

	//In:
	//	3 x 1 normal vector to the circle-plane
	//	3 x 1 direction of interest u
	//	3 x 3 output rotation matrix for that direction


	double add[3];
	double prod[3][3];

	add[0] = z[0] + u[0];
	add[1] = z[1] + u[1];
	add[2] = z[2] + u[2];

	matMulPoints(add, add, prod);
	double div = mulPoints(z, u, 3) + 1;

	//printf("Denom: %f\n", div);

	for(int i = 0; i < 3; ++i){
		for(int j = 0; j < 3; ++j){
			prod[i][j] = prod[i][j] / div;
		}
	}
	prod[0][0] --; //Subtract I from the matrix
	prod[1][1] --;
	prod[2][2] --;

	for(int i = 0; i < 3; ++i) {
		for(int j = 0; j < 3; ++j){
			out[i][j] = prod[i][j];
		}
	}

}

double getKernel(double* n1, double* n2){
	//printf("getKernel\n");
	double width = 1; //WIDTH PARAMETER, LEFT CONSTANT FOR NOW

	double mul = mulPoints(n1,n2,3);

	if(mul < 0){
		mul *= -1;
	}

	if(mul > 1 && mul <= 1.0001){ //Catch any floating point errors
		mul = 1;
	}
	else if(mul < -1 && mul >= -1.0001){
		mul = -1;
	}

	double d = acos(mul);


	//spherical Gaussian function
	
	double nd2 = (d*d);
	double out = exp(nd2/(width*width));

	//printf("\t\tkernel: %f\n",out);
	return out;
}