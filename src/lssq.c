#include "lssq.h"

LSSQ_ARENA *lssq_initialize(SYMMAT33 **b_matrices, int nmatrices, int nentries)
{
    int i;
    int ncol = nentries;
    int nrow = nmatrices;
    int r,c;
    int gj_retval;
    
    LSSQ_ARENA *lssq = malloc(sizeof(LSSQ_ARENA));
    
    SYMMAT33 **bmat_tmp = malloc(sizeof(SYMMAT33 *) * nrow);
    for (i=0; i<nrow; i++)  {
        bmat_tmp[i] = malloc(sizeof(SYMMAT33));
        for (c=0; c<ncol; c++) {
            bmat_tmp[i]->data[c] = -1.0e+3 * b_matrices[i]->data[c];
        }
    }
    
    double *xmean   = malloc(sizeof(double)*ncol);
    memset(xmean, 0, ncol*sizeof(double));
    
    double *sigmax = malloc(sizeof(double) * ncol);
    memset(sigmax, 0, ncol*sizeof(double));
    
    double **btb_inv = malloc(sizeof(double *) * ncol);
    double **btb = malloc(sizeof(double *) * ncol);
    double **btb_copy = malloc(sizeof(double *) * ncol);
    for (i=0; i<ncol; i++) {
        btb[i] = malloc(sizeof(double) * ncol);
        btb_inv[i] = malloc(sizeof(double) * ncol);
        btb_copy[i] = malloc(sizeof(double) * ncol);
    }
    
    for (r=0; r<nrow; r++) {
        for (c=0; c<ncol; c++) {
            xmean[c] += (bmat_tmp[r]->data[c])/((double)nrow);
        }
    }
    
    for (r=0; r<nrow; r++) {
        for (c=0; c<ncol; c++) {
            bmat_tmp[r]->data[c] -= xmean[c];
        }
    }
    
    for (r=0; r<nrow; r++) {
        for (c=0; c<ncol; c++) {
            sigmax[c] += (bmat_tmp[r]->data[c])*(bmat_tmp[r]->data[c]);
        }
    }
    
    for (c=0; c<ncol; c++)
        sigmax[c] = sqrt(sigmax[c]/((double)(nrow-1)));
    
    for (r=0; r<ncol; r++) {
        
        for (c=0; c<ncol; c++) {
            
            int i = 0;
            
            if (r == c) {
                btb_inv[r][c] = 1;
            } else {
                btb_inv[r][c] = 0;
            }
            
            for (i=0; i<nrow; i++) {
                
                double tmp1 = 0;
                double tmp2 = 0;
                tmp1 = (bmat_tmp[i]->data[r]);
                tmp2 = (bmat_tmp[i]->data[c]);
                
                btb[r][c] += tmp1 * tmp2;
                
            }
            
        }
        
    }
    
    for (r=0; r<ncol; r++) {
        
        for (c=0; c<ncol; c++) {
            
            btb[r][c] = btb[r][c] / ( ((double)(nrow-1)) * sigmax[r] * sigmax[c] );
            /* printf("%e  ", ( ((double)(nrow-1)) * sigmax[r] * sigmax[c] )); */
            
        }
        
        /* printf("\n"); */
        
    }
    
    for (r=0; r<ncol; r++) {
        for (c=0; c<ncol; c++) {
            btb_copy[r][c] = btb[r][c];
        }
    }
    
    gj_retval = lssq_gauss_jordan_solve(btb_copy, btb_inv, NULL, NULL, ncol);
    
    /* retval will be checked shortly */
    lssq->bmat_tmp     = bmat_tmp;
    lssq->xmean        = xmean;
    lssq->sigmax       = sigmax;
    lssq->ncols        = ncol;
    lssq->nrows        = nrow;
    lssq->btranspose_b = btb;
    lssq->btranspose_b_inv = btb_inv;
    
    if (gj_retval == 0) {
        fprintf(stderr, "lssq_initialize: System is not invertable.\n");
        lssq_free(lssq);
        lssq = NULL;
        /* so that null is returned */
    }
	
    for (r=0; r<ncol; r++) {
        free(btb_copy[r]);
    }
    free(btb_copy);
    
    return lssq;
    
    /*************************************************************************/
    
}

void lssq_free(LSSQ_ARENA *lssq)
{
    int i;
    
    for (i=0; i<lssq->nrows; i++)
        free(lssq->bmat_tmp[i]);
    
    free(lssq->bmat_tmp);
    
    free(lssq->xmean);
    free(lssq->sigmax);
    
    for (i=0; i<lssq->ncols; i++) {
        free(lssq->btranspose_b[i]);
        free(lssq->btranspose_b_inv[i]);
    }
    
    free(lssq->btranspose_b);
    free(lssq->btranspose_b_inv);
    free(lssq);
    
}

double lssq_compute(LSSQ_ARENA *lssq, double *X, double *B)
{
    
    int r = 0;
    int c = 0;
    
    SYMMAT33 **bmat_tmp = lssq->bmat_tmp;
    double nobs         = (double)(lssq->nrows);
    double *sigmax      = lssq->sigmax;
    double *xmean       = lssq->xmean;
    double ymean        = 0;
    double sigmay       = 0;
    
    for (r=0; r<lssq->nrows; r++) {
        B[r] = log(B[r]);
        ymean += B[r]/nobs;
    }
    
    for (r=0; r<lssq->nrows; r++) {
        B[r] -= ymean;
        sigmay += B[r]*B[r]/(nobs-1);
    }
    
    sigmay = sqrt(sigmay);
    
    /* b_transpose times B */
    for (r=0; r<6; r++) {
        X[r] = 0;
        for (c=0; c<lssq->nrows; c++) {
            X[r] += B[c]*(bmat_tmp[c]->data[r]);
        }
        X[r] /= ( sigmay * sigmax[r] * (nobs-1) );
    }
    
    /* btb_inverse times (b_transpose*B) */
    /* since we are finished with the original data, I am using B to store
     temporarily */
    for (r=0; r<6; r++) {
        double tmp = 0;
        for (c=0; c<6; c++) {
            tmp += X[c]*lssq->btranspose_b_inv[c][r];
        }
        B[r] = tmp;
    }
    
    /* move data back to X */
    /* for (r=0; r<6; r++) X[r] = B[r]; */
    
    double S0 = 0;
    
    for (r=0; r<6; r++) {
        X[r] = B[r] * (sigmay/sigmax[r]);
        S0 += X[r] * xmean[r];
    }
    
    /*
     for (r=0; r<6; r++) {
     S0 += X[r] * xmean[r];
     }
     */
    
    S0 = exp(ymean - S0);
    
    return S0;
    
}

int lssq_gauss_jordan_solve (double **A, double **I, double *x, double *b, int n)
{
    
    int pivot_row = 0;
    int pivot_col = 0;
    int r,c;
    double *piv_tmp;
    
    for (r=pivot_row; r<n; r++) {
        
        double div = A[r][pivot_col];
        int r1;
        
        if (div == 0) {
            /* perform row swap with first row with nonzero leading coef. */
            int r2;
            
            for (r2=r+1; r2<n; r2++) {
                if (A[r2][pivot_col] != 0) {
                    
                    piv_tmp = A[r];
                    A[r]    = A[r2];
                    A[r2]   = piv_tmp;
                    
                    piv_tmp = I[r];
                    I[r]    = I[r2];
                    I[r2]   = piv_tmp;
                    
                    div = A[r][pivot_col];
                    
                    break;
                }
            }
            
            if (r2 == n) {
                /* could not find suitable row. */
                fprintf(stderr, "lssq_gauss_jordan_solve: Singular Matrix Detected.\n");
                return 0;
            }
        }
        
        /* divide pivot row by leading coef of pivot col */
        for (c=pivot_col; c<n; c++) {
            A[r][c] /= div;
        }
        for (c=0; c<n; c++) {
            I[r][c] /= div;
        }
        //b[r] /= div;
        
        /* now perform row operations to eliminate the other coefs */
        for (r1=pivot_row; r1<n; r1++) {
            
            /* don't mess with the pivot row */
            if (r1 == r) continue;
            
            if (A[r1][pivot_col] == 0) continue;
            
            double coef = -(A[r1][pivot_col]);
            
            
            for (c=pivot_col; c<n; c++) {
                A[r1][c] += coef*A[r][c];
            }
            for (c=0; c<n; c++) {
                I[r1][c] += coef*I[r][c];
            }
            
            //b[r1] += coef*b[r];
        }
        
        pivot_col++;
        
    }
    
    return 1;
    
}

/*
 Computes all eigenvalues and eigenvectors of a real symmetric matrix a[1..n][1..n]. On
 output, elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a.
 v[1..n][1..n] is a matrix whose columns contain, on output, the normalized eigenvectors of
 a. nrot returns the number of Jacobi rotations that were required.
 */

#ifdef LSSQ_ROTATE
#undef LSSQ_ROTATE
#endif

#define LSSQ_ROTATE(a,i,j,k,l) g = a[i][j];             \
h = a[k][l];             \
a[i][j] = g-s*(h+g*tau); \
a[k][l] = h+s*(g-h*tau);

void lssq_jacobi(float **a, int n, float d[], float **v, int *nrot,
                 float *b, float *z)
{
    int j,iq,ip,i;
    float tresh,theta,tau,t,sm,s,h,g,c;//,*b,*z;
    
    /*
     b = lssq_vector(1,n);
     z = lssq_vector(1,n);
     */
    for (ip=1;ip<=n;ip++) { 
        /* Initialize to the identity matrix. */
        for (iq=1;iq<=n;iq++) {
            v[ip][iq] = 0.0;
        }
        v[ip][ip] = 1.0;
    }
    
    for (ip=1;ip<=n;ip++) { 
        /* Initialize b and d to the diagonal of a. */
        b[ip] = d[ip] = a[ip][ip]; 
        z[ip] = 0.0; /* This vector will accumulate terms of
                      the form tapq as in equation (11.1.14). */
    }
    
    *nrot = 0;
    for (i=1;i<=50;i++) {
        
        sm = 0.0;
        
        for (ip=1;ip<=n-1;ip++) { 
            /* Sum o-diagonal elements. */
            for (iq=ip+1;iq<=n;iq++) {
                sm += fabs(a[ip][iq]);
            }
        }
        
        if (sm == 0.0) { 
            /* The normal return, which relies on quadratic convergence to machine underflow. */
            /*
             lssq_free_vector(z,1,n);
             lssq_free_vector(b,1,n);
             */
            lssq_eigen_sort(d, v, n);
            return;
        }
        
        if (i < 4) {
            tresh = 0.2*sm/(n*n); /* ...on the first three sweeps. */
        } else {
            tresh = 0.0;          /* ...thereafter.                */
        }
        
        for (ip=1;ip<=n-1;ip++) {
            
            for (iq=ip+1;iq<=n;iq++) {
                
                g = 100.0*fabs(a[ip][iq]);
                
                /* After four sweeps, skip the rotation if the o-diagonal element is small. */
                
                if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
                    && (float)(fabs(d[iq])+g) == (float)fabs(d[iq])) {
                    a[ip][iq] = 0.0;
                } else if (fabs(a[ip][iq]) > tresh) {
                    
                    h = d[iq]-d[ip];
                    
                    if ((float)(fabs(h)+g) == (float)fabs(h)) {
                        t = (a[ip][iq])/h;
                    } else {
                        theta = 0.5*h/(a[ip][iq]); /*  Equation (11.1.10). */
                        t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                        if (theta < 0.0) 
                            t = -t;
                    }
                    
                    c   = 1.0/sqrt(1+t*t);
                    s   = t*c;
                    tau = s/(1.0+c);
                    h   = t*a[ip][iq];
                    
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq] = 0.0;
                    
                    for (j=1;j<=ip-1;j++) { 
                        /* Case of rotations 1 <= j < p.*/
                        LSSQ_ROTATE(a,j,ip,j,iq)
                    }
                    
                    for (j=ip+1;j<=iq-1;j++) { 
                        /* Case of rotations p < j < q. */
                        LSSQ_ROTATE(a,ip,j,j,iq)
                    }
                    
                    for (j=iq+1;j<=n;j++) { 
                        /* Case of rotations q < j <= n. */
                        LSSQ_ROTATE(a,ip,j,iq,j)
                    }
                    
                    for (j=1;j<=n;j++) {
                        LSSQ_ROTATE(v,j,ip,j,iq)
                    }
                    
                    ++(*nrot);
                }
            }
        }
        
        for (ip=1;ip<=n;ip++) {
            b[ip] += z[ip];
            d[ip] = b[ip];   /* Update d with the sum of tapq, */
            z[ip] = 0.0;     /* and reinitialize z. */
        }
        
    }
    
    fprintf(stderr, "Too many iterations in routine jacobi\n");
}

void lssq_eigen_sort(float d[], float **v, int n) 
/* Given the eigenvalues d[1..n] and eigenvectors v[1..n][1..n] as output
 from jacobi (ss11.1) or tqli (ss11.3), this routine sorts the eigenvalues
 into descending order, and rearranges the columns of v correspondingly.
 The method is straight insertion. */
{
    int k,j,i;
    float p;
    
    for (i=1;i<n;i++) {
        
        p = d[k = i];
        
        for (j=i+1;j<=n;j++)
            if (d[j] >= p) p = d[k = j];
        
        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            
            for (j=1;j<=n;j++) {
                p = v[j][i];
                v[j][i] = v[j][k];
                v[j][k] = p;
            }
        }
    }
}

float *lssq_vector(int nl, int nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
    float *v;
    v=(float *)malloc((size_t) ((nh-nl+1+1)*sizeof(float)));
    if (!v) fprintf(stderr, "allocation failure in vector()");
    return v-nl+1;
}

void lssq_free_vector(float *v, int nl, int nh)
/* free a float vector allocated with vector() */
{
    free((char *) (v+nl-1));
}
