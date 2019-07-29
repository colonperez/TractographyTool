/*
 *  mow_recon.h
 *  track_tools
 *
 *  Created by William Triplett on 6/15/10.
 *
 */

#include "nifti1_io.h"
#include "track_track.h"
#include "matrices.h"
#include "lssq.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _WIN32
const char DIRSEP = '\\';
#else
const char DIRSEP = '/';
#endif

#define DO_TENSOR 0

const double DECO_P        = 2.0;
const double DECO_EVALS[3] = {1.5, 0.4, 0.4};
const double DELTA_SM      = 16.9100;
const double DELTA_LG      = 29.7300;
const double DIFF_RADIUS   = 10.0;
const double PROB_THRESH   = 0.5;

/****************************************************************************
 * Data Structures
 ****************************************************************************/

typedef struct DIFF_DATA {
    nifti_brick_list *nii_brick_list;
    nifti_image *nii_image;
    nifti_image *mask;
    nifti_image *S0;
    double *single_voxel_storage; 
    int n_volumes;
    float *bvals;
    float *bvecs;
    SYMMAT33 **b_matrices;
    int *b_high_ind;
    int *b_low_ind;
    int n_b_high;
    int n_b_low;
    double delta_lg;
    double delta_sm;
} DIFF_DATA;    

typedef struct ICOS_TESS {
    int num_vertices;
    float **vertices;
    int **connectivity;
} ICOS_TESS;

typedef struct MOW_RECON {
    DIFF_DATA *diff;
    ICOS_TESS *reco_tess;
    ICOS_TESS *deco_tess;
    ICOS_TESS *restart_tess;
    double deco_p;
    double deco_evals[3];
    double diff_radius;
    double *A_matrix;
    double *reco_matrix;
    double **Sig_i;
    double prob_thresh;
    SYMMAT33 **D_i_inv;
    char *output_directory;
    char *data_filename;
    char *mask_filename;
    char *bval_filename;
    char *bvec_filename;
    char *S0_filename;
    char *datadir;
    int S0compute;
    int log_bad_voxels;
    int num_output_files;
} MOW_RECON;

typedef struct OUTPUT_DATA {
    nifti_image **nim;
    nifti_brick_list **nbl;
    int num_images;
} OUTPUT_DATA;

typedef struct MAXIMA {
    double value;
    int index;
} MAXIMA;

/****************************************************************************
 * File I/O routines
 ****************************************************************************/

ICOS_TESS *load_tess_from_file(const char *tess_file);
void read_bvals_from_file(char *filename, DIFF_DATA *diff);
void read_bvecs_from_file(char *filename, DIFF_DATA *diff);
void read_diff_data_from_file(char *filename, DIFF_DATA *diff);
int read_acsii_file_to_float_array(char *filename, float *data, int data_size);
int save_output (const char *basedir, OUTPUT_DATA *output);

/****************************************************************************
 * Optimization / Maxima Finding routines
 ****************************************************************************/

int find_local_maxima(ICOS_TESS *domain, double *values, double min_values,
		      ICOS_TESS *restart, MAXIMA *maxima_list);    
int maxima_compare(const void *m1, const void *m2);

/****************************************************************************
 * Core processing routines
 ****************************************************************************/

void mow_initialize_opts(MOW_RECON *mow, int argc, char **argv);
void mow_print_usage(void);
float *compute_S0 (DIFF_DATA *diff, LSSQ_ARENA *lssq, int log_bad_voxels);
float compute_FA_from_evals(float evals[]);
void compute_bmatrix(DIFF_DATA *diff);
double *make_A_matrix(MOW_RECON *mow);
double *make_recon_matrix(MOW_RECON *mow);
void make_reconstruction_matrix(DIFF_DATA *diff, ICOS_TESS *tess);
double *compute_reconstruction_weights(MOW_RECON *mow, int x, int y, int z);

/****************************************************************************
 * Data handling
 ****************************************************************************/

OUTPUT_DATA *initialize_output (nifti_image *template, int nbricks);
void add_maxima_to_output(OUTPUT_DATA *output, int x, int y, int z,
			  float **vertlist, MAXIMA *maxima_list, int n_maxima);
void replace_realsmall_w_zeros(double *input, int size, double tol);
void __tokenize_line(char *line, const char *tokens, float *split, int *nsplit);
int load_voxel_double_highb(DIFF_DATA *diff, int x, int y, int z);
int load_voxel_double_all(DIFF_DATA *diff, int x, int y, int z, double *dest);
double read_nii_voxel_anytype(void *src, int index, int datatype);



