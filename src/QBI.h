#ifndef QBI_H_
#define QBI_H_

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdbool.h>

#include"qbi_nifti1_io.h"
#include"qbi_ConvolutionAndPrediction.h"
#include"qbi_Matrices.h"


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

typedef struct ICOS_TESS {
    int num_vertices;
    float **vertices;
    int **connectivity;
} ICOS_TESS;


typedef struct OUTPUT_DATA {
    nifti_image **nim;
    nifti_brick_list **nbl;
    int num_images;
} OUTPUT_DATA;

typedef struct MAXIMA {
    double value;
    int index;
} MAXIMA;

typedef struct SYMMAT33 {
    float data[6];
} SYMMAT33;

typedef struct MAT33 {
    float data[9];
} MAT33;

typedef struct DIFF_DATA {
    nifti_brick_list *nii_brick_list;
    nifti_image *nii_image;
    nifti_image *mask;
    nifti_image *S0;
    nifti_image *GFA;
    double *single_voxel_storage;
    int n_volumes;
    float *bvecs;
    float *bvals;
    int *b_high_ind;
    int *b_low_ind;
    int n_b_high;
    int n_b_low;

} DIFF_DATA;

typedef struct QBI_RECON {
    DIFF_DATA *diff;
    ICOS_TESS *reco_tess;
    ICOS_TESS *restart_tess;
    double deco_p;
    double deco_evals[3];
    double prob_thresh;
    char *output_directory;
    char *data_filename; //e
    char *mask_filename;
    char *bvec_filename; //Q
    char *bval_filename;
    char *GFA_filename;
    int GFAcompute;
    char *datadir;
    int log_bad_voxels;
    int num_output_files;
    char *S0_filename;
} QBI_RECON;

void qbi_initialize_opts(QBI_RECON *qbi, int argc, char **argv); //Reads in command line arguments and initializes all of the structs and variables
ICOS_TESS *load_tess_from_file(const char *tess_file); //Loads in spherical tessellation file
void read_diff_data_from_file(char *filename, DIFF_DATA *diff); // Reads in diffusion data
void __tokenize_line(char *line, const char *tokens, float *split, int *nsplit); //Used when reading in text for 
int read_acsii_file_to_float_array(char *filename, float *data, int data_size); //Reads in a file containing numbers and creates a list of floats with the read-in data
int maxima_compare(const void *m1, const void *m2); //Used by MOW to compare maximas, checks for duplicate maximas, not used in my code
int find_local_maxima(ICOS_TESS *domain, double *values, double min_value, ICOS_TESS *restart, MAXIMA *maxima_list); //Used by MOW to find maximas, relies on the tess files to find maxima, not used in my code
void qbi_print_usage(void); //Prints the arguments used in the program, I edited this to only use arguments used in QBI
OUTPUT_DATA *initialize_output (nifti_image *template, int nfiles); //Initializes the output data structure
int nii_voxel3_index(nifti_image *nim, const int x, const int y, const int z); //Used in converting nifti files to the data structs
void read_bvecs_from_file(char *filename, DIFF_DATA *diff); //Reads in the bvecs from the bvec file into a single one dimensional array
int load_voxel_double_all(DIFF_DATA *diff, int x, int y, int z, double *dest); //Loads a single voxel from the data
int load_voxel_double_highb(DIFF_DATA *diff, int x, int y, int z); //Loads a single voxel from the data and divides it by that voxel's S0
void add_maxima_to_output(OUTPUT_DATA *output, int x, int y, int z, double **vertlist, MAXIMA *maxima_list, int n_maxima); //Adds the list of maxima found into the output data struct
double read_nii_voxel_anytype(void *src, int index, int datatype); //Used in converting nifti files to the data structs
int nii_recast_to_int32 (nifti_image *nim); //Used in converting nifti files to the data structs
int save_output (const char *basedir, OUTPUT_DATA *output); //Saves output data structure to files
void read_bvals_from_file(char *filename, DIFF_DATA *diff);

void getA(double ** Q, double ** U, double ** V, int k, int m, int n, int p, double** A);
void QBall(double* e, double** A, int n, int m, double* ODF);



int grad_find_max(double* ODF, int n, MAXIMA* maxima_list, int adj[50][4]);
int temp_find_max(double* ODF, int ODFLen, MAXIMA* maxima_list);

#endif
