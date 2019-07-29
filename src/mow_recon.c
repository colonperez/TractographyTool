/*
 *  mow_recon.c
 *  track_tools
 *
 *  Created by William Triplett on 6/15/10.
 *
 */

#include "../../TrackTools_GUI/tt_current_version.h"
#include "mow_recon.h"
#include "nnls.h"
#include "lssq.h"

int main (int argc, char **argv)
{
    
    MOW_RECON mow;
    DIFF_DATA diff;
    
    mow.diff = &diff;
    
    mow_initialize_opts(&mow, argc, argv);
    
    /* Load the precomputed spherical domains */
    char *strbuf = malloc(sizeof(char)*strlen(mow.datadir) + 15);
    
    fprintf(stderr, "mow_recon (%d) starting...\n", kTT_CURRENT_VERSION);
    
    fprintf(stderr, "Loading spherical domains...\n");
    sprintf(strbuf, "%s%c%s", mow.datadir, DIRSEP, "tess_L1.dat");
    ICOS_TESS *restart_tess = load_tess_from_file(strbuf);
    
    sprintf(strbuf, "%s%c%s", mow.datadir, DIRSEP, "tess_L2.dat");
    ICOS_TESS *deco_tess    = load_tess_from_file(strbuf);
    
    sprintf(strbuf, "%s%c%s", mow.datadir, DIRSEP, "tess_L3.dat");
    ICOS_TESS *reco_tess    = load_tess_from_file(strbuf);
    if (reco_tess == NULL || deco_tess == NULL || restart_tess == NULL) {
        fprintf(stderr, "Spherical tessellation files could not be loaded.\n");
        fprintf(stderr, "Make sure that the -datadir option points to the\n");
        fprintf(stderr, "directory that contains the 'tess_L*.dat' files.\n");
        exit(1);
    }
    free(strbuf); strbuf = NULL;
    
    mow.reco_tess = reco_tess;
    mow.deco_tess = deco_tess;    
    mow.restart_tess = restart_tess;
    
    /* Set up the output data structure that will become the NIFTI direction 
     files. */
    fprintf(stderr, "Initializing output data structures...\n");
    OUTPUT_DATA *output = initialize_output(diff.nii_image, mow.num_output_files);
    
    /* self-explanatory */
    fprintf(stderr, "Setting up B matrix...\n");
    compute_bmatrix(&diff);
    fprintf(stderr, "Computing A matrix (deconvolution)...\n");
    make_A_matrix(&mow);
    fprintf(stderr, "Computing reconstruction matrix...\n");
    make_recon_matrix(&mow);
    
    int n_reco_dirs = reco_tess->num_vertices;
    int n_deco_dirs = deco_tess->num_vertices;
    
    /* set up the NNLS solver for multiple solutions using the same A matrix. */
    fprintf(stderr, "Preparing non-negative least squares solver...\n");
    NNLS_ARENA *nnls = nnls_initialize(mow.A_matrix, diff.n_b_high, n_deco_dirs);
    
    /* initialize the coefficient storage for a single voxel */
    double *coef = malloc(sizeof(double) * n_reco_dirs);
    if (coef == NULL) {
        fprintf(stderr, "Unable to allocate memory for coefficients.\n");
        exit(1);
    }
    
    /* set up the maxima list. */
    MAXIMA *maxima_list  = malloc(sizeof(MAXIMA) * n_reco_dirs);
    if (maxima_list == NULL) {
        fprintf(stderr, "Unable to allocate memory for maxima list.\n");
        exit(1);
    }
    
    /*************************************************************************/
    
    if (mow.S0compute == 1) {
        
        LSSQ_ARENA *lssq = lssq_initialize(diff.b_matrices, diff.n_volumes, 6);
        if (lssq == NULL) {
            fprintf(stderr, "Unable to initialize least-squares solver.\n");
            exit(1);
        }
        
        float *S0_image = compute_S0(&diff, lssq, mow.log_bad_voxels);
        nifti_image *S0_nim = nifti_simple_init_nim();
        memcpy(S0_nim, diff.nii_image, sizeof(nifti_image));
        
        S0_nim->datatype = DT_FLOAT32;
        S0_nim->ndim     = 3;
        S0_nim->nbyper   = 4;
        S0_nim->nt       = 1;
        S0_nim->nvox     = S0_nim->nx * S0_nim->ny * S0_nim->nz;
        S0_nim->dim[4]   = 1;
        S0_nim->dim[0]   = 3;
        S0_nim->data     = S0_image;
        S0_nim->fname    = mow.S0_filename;
        S0_nim->iname    = mow.S0_filename;
        S0_nim->cal_max  = 0.0;
        S0_nim->cal_min  = 0.0;
        sprintf(S0_nim->descrip, "TrackTools MOW S0 Data");
        diff.S0 = S0_nim;
        lssq_free(lssq);
        lssq = NULL;
        
        if (mow.S0_filename != NULL) {
            fprintf(stderr, "Saving S0 image to %s.\n", mow.S0_filename);
            znzFile fp = znzopen(mow.S0_filename, "wb", 0);
            nifti_image_write_hdr_img2(S0_nim, 1, "wb", fp, NULL);
        }
    }
    
    fprintf(stderr, "Starting MOW reconstruction...\n");
    /*************************************************************************/
    
    int vx, vy, vz;    
    double diff_time = diff.delta_lg - diff.delta_sm/3.0;
    double determ_d  = (mow.deco_evals[0] *
                        mow.deco_evals[1] *
                        mow.deco_evals[2]);
    double w_scale = sqrt( pow((4.0*M_PI*diff_time),3) * determ_d );
    double *reco_matrix = mow.reco_matrix;
    
    for (vz=0; vz<diff.nii_image->nz; vz++) {
        for (vy=0; vy<diff.nii_image->ny; vy++) {
            for (vx=0; vx<diff.nii_image->nx; vx++) {
                
                double min = 1.0e+99, max = 0;
                double *x = NULL; /* nnls solution vector */
                double *data = NULL; /* diff-weighted data */
                int dec;
                int rec;
                int n_maxima = 0;
                
                int load_ok = load_voxel_double_highb(&diff, vx, vy, vz);
                
                if (-1 == load_ok) {
                    if (mow.log_bad_voxels != 0) {
                        fprintf(stderr, "  WARNING: Voxel [%d,%d,%d] was not reconstructed (S0=0).\n",
                                vx, vy, vz);
                    }
                    continue;
                } else if (-2 == load_ok) {
                    if (mow.log_bad_voxels != 0) {
                        fprintf(stderr, "  WARNING: Voxel [%d,%d,%d] was not reconstructed (nan/inf).\n",
                                vx, vy, vz);
                    }
                    continue;
                } else if (0 == load_ok) {
                    continue;
                }
                
                data = diff.single_voxel_storage;
                nnls_compute(nnls, mow.A_matrix, data, 1);
                x = nnls->x;
                
                /* reset coef to 0 for next run */
                memset(coef, 0, n_reco_dirs*sizeof(double));
                
                /* reconstruct coefficients for vertices */
                for (rec=0; rec<n_reco_dirs; rec++) {
                    
                    /* Eqn #16 */
                    for (dec=0; dec<n_deco_dirs; dec++) 
                        coef[rec] += (x[dec]/w_scale) * reco_matrix[dec*n_reco_dirs+rec];
                    
                    if (coef[rec] > max) 
                        max = coef[rec];
                    if (coef[rec] < min)
                        min = coef[rec];
                    
                }
                
                /* does min/max scaling. Probably not necessary for maxima finding. */
                if (1) {
                    for (rec=0; rec<n_reco_dirs && min != max; rec++)
                        coef[rec] = (coef[rec]-min)/(max-min);
                }
                
                n_maxima = find_local_maxima(reco_tess, coef, mow.prob_thresh,
                                             restart_tess, maxima_list);
                
                add_maxima_to_output(output, vx, vy, vz,
                                     reco_tess->vertices, maxima_list, n_maxima);
                
            }
        }
        
        fprintf(stderr, "Slice: %d of %d Complete.\n", vz, diff.nii_image->nz);
        fflush(stderr);
    }
    
    fprintf(stderr, "MOW Reconstruction complete... saving output...\n");
    save_output(mow.output_directory, output);
    
    fprintf(stderr, "Done.\n");
    return 0;
}

/*****************************************************************************
 * computes the FA given 3 eigenvalues.
 *****************************************************************************/

float compute_FA_from_evals(float evals[])
{
    float mean_ev;
    float dev[3];
    float denom = 0.0;
    float FA = 0.0;
    int i;
    
    mean_ev = (evals[0]+evals[1]+evals[2])/3.0;
    
    for (i=0; i<2; i++) {
        dev[i] = (evals[i]-mean_ev);
        dev[i] *= dev[i];
        FA += dev[i];
        denom += (evals[i]*evals[i]);
    }
    
    FA = sqrt( (3.0 * FA)/(2.0 * denom) );
    
    return FA;
    /* return (FA > 1.0) ? 1.0 : FA; */
    
}

/*****************************************************************************
 * computes an S0 image using DTI method.
 *****************************************************************************/

float *compute_S0 (DIFF_DATA *diff, LSSQ_ARENA *lssq, int log_bad_voxels)
{
    
    int vx=0;
    int vy=0;
    int vz=0;
    int nx = diff->nii_image->nx;
    int ny = diff->nii_image->ny;
    int nz = diff->nii_image->nz;
    
    double *X = malloc(sizeof(double) * 6);
    float *S0_image = malloc(sizeof(float)*nx*ny*nz);
    if (X == NULL || S0_image == NULL) {
        fprintf(stderr, "Unable to allocate memory for S0 image.\n");
        return NULL;
    } else {
        fprintf(stderr, "Computing S0 image...\n");
    }
    
    memset(S0_image, 0, sizeof(float)*nx*ny*nz);
    
#ifdef DO_TENSOR    
    float **tensor = malloc(sizeof(float *) * 4);
    float **evecs  = malloc(sizeof(float *) * 4);
    float *b       = malloc(sizeof(float)*4);
    float *z       = malloc(sizeof(float)*4);
    for (vz=0; vz<4; vz++) {
        tensor[vz] = malloc(sizeof(float)*4);
        evecs[vz]  = malloc(sizeof(float)*4);
        memset(tensor[vz], 0, sizeof(float)*4);
        memset(evecs[vz], 0, sizeof(float)*4);
    }
#endif
    
    for (vz=0; vz<nz; vz++) {
        for (vy=0; vy<ny; vy++) {
            for (vx=0; vx<nx; vx++) {
                
                double S0;
#ifdef DO_TENSOR		
                int nrot = 0;
                float d[4];
                float FA;
#endif		 
                /*
                 if (vx == 54 && vy == 67 && vz == 35) {
                 printf(".\n");
                 }
                 */
                
                int index = nii_voxel3_index(diff->nii_image, vx, vy, vz);
                
                int load_ok = load_voxel_double_all(diff, vx, vy, vz, NULL);
                
                if (-1 == load_ok) {
                    if (log_bad_voxels != 0) {
                        fprintf(stderr, "  WARNING: Voxel [%d,%d,%d] was not loaded.\n",
                                vx, vy, vz);
                    }
                    continue;
                } else if (0 == load_ok) {
                    /* mask hit */
                    continue;
                }
                
                S0 = lssq_compute(lssq, X, diff->single_voxel_storage);
                
                if (isnan(S0) || isinf(S0)) {
                    if (log_bad_voxels != 0) {
                        fprintf(stderr, "  WARNING: Voxel [%d,%d,%d] nan/inf in S0.\n",
                                vx, vy, vz);
                    }
                    continue;
                }
                
                S0_image[index] = (float) S0;
                
#ifdef DO_TENSOR
                tensor[1][1] = (float)X[0];
                tensor[2][2] = (float)X[1];
                tensor[3][3] = (float)X[2];
                tensor[1][2] = tensor[2][1] = (float)X[3];
                tensor[1][3] = tensor[3][1] = (float)X[4];
                tensor[2][3] = tensor[3][2] = (float)X[5];
                
                lssq_jacobi(tensor, 3, d, evecs, &nrot, b, z);
                FA = compute_FA_from_evals(d);
                
                
                /* evals[0 ..3] = d[1 ... 3], evecs[1..3][1..3] */
                tensor[1][1] = tensor[1][1];
#endif
                
            }
        }
    }
    
#ifdef DO_TENSOR
    free(b);
    free(z);
#endif
    free(X);
    
    return S0_image;
    
}

/*****************************************************************************
 * reads command line options and sets up the preferences accordingly. 
 *****************************************************************************/

void mow_initialize_opts(MOW_RECON *mow, int argc, char **argv)
{
    int opt;
    
    /*
     char *direction_files_multi[5] = {
     "V_00_all.nii.gz", 
     "V_01_all.nii.gz",
     "V_02_all.nii.gz",
     "V_03_all.nii.gz",
     "V_04_all.nii.gz"};
     */
    
    /* these may be overwritten by command line options */
    mow->deco_p           = (double) DECO_P;
    mow->diff_radius      = DIFF_RADIUS;
    mow->deco_evals[0]    = (double) DECO_EVALS[0];
    mow->deco_evals[1]    = (double) DECO_EVALS[1];
    mow->deco_evals[2]    = (double) DECO_EVALS[2];
    mow->prob_thresh      = (double) PROB_THRESH;
    mow->num_output_files = 5;
    mow->S0compute        = 0;
    mow->log_bad_voxels   = 0;
    mow->datadir          = ".";
    mow->diff->delta_lg   = (double) DELTA_LG;
    mow->diff->delta_sm   = (double) DELTA_SM;
    mow->data_filename    = NULL;
    mow->bval_filename    = NULL;
    mow->bvec_filename    = NULL;
    mow->mask_filename    = NULL;
    mow->S0_filename      = NULL;
    
    for (opt=1; opt<argc; opt++) {
        if (0 == strcmp(argv[opt], "-h")) {
            mow_print_usage();
            exit(0);
        }
        if (0 == strcmp(argv[opt], "-odir")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -odir requires an argument.\n");
                exit(1);
            }
            mow->output_directory = argv[opt+1];
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-S0")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -S0 requires an argument.\n");
                exit(1);
            }
            mow->S0_filename = argv[opt+1];
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-S0compute")) {
            mow->S0compute = 1;
            continue;
        } else if (0 == strcmp(argv[opt], "-log-bad-voxels")) {
            mow->log_bad_voxels = 1;
            continue;
        } else if (0 == strcmp(argv[opt], "-pthresh")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -pthresh requires an argument.\n");
                exit(1);
            }	    
            mow->prob_thresh = fabs(atof(argv[opt+1]));
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-mask")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -mask requires an argument.\n");
                exit(1);
            }	    
            mow->mask_filename = argv[opt+1];
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-datadir")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -datadir requires an argument.\n");
                exit(1);
            }	    
            mow->datadir = argv[opt+1];
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-ndir")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -ndir requires an argument.\n");
                exit(1);
            }	    
            mow->num_output_files = atoi(argv[opt+1]);
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-data")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -data requires an argument.\n");
                exit(1);
            }	    
            mow->data_filename = argv[opt+1];
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-bval")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -bval requires an argument.\n");
                exit(1);
            }	    
            mow->bval_filename = argv[opt+1];
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-bvec")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -bvec requires an argument.\n");
                exit(1);
            }	    
            mow->bvec_filename = argv[opt+1];
            opt++;
            continue;	    
        } else if (0 == strcmp(argv[opt], "-radius")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -radius requires an argument.\n");
                exit(1);
            }	    
            mow->diff_radius = atof(argv[opt+1]);
            opt++;
            continue;	    
        } else if (0 == strcmp(argv[opt], "-delta_lg")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -delta_lg requires an argument.\n");
                exit(1);
            }	    
            mow->diff->delta_lg = atof(argv[opt+1]);
            opt++;
            continue;	    
        } else if (0 == strcmp(argv[opt], "-delta_sm")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -delta_sm requires an argument.\n");
                exit(1);
            }	    
            mow->diff->delta_sm = atof(argv[opt+1]);
            opt++;
            continue;	    
        } else {
            fprintf(stderr, "Ignoring junk on command line: %s\n", argv[opt]);
            fprintf(stderr, "Use -h to see a list of options.\n");
            exit(1);
        }
    }
    
    if (mow->S0compute == 1) {
        if (mow->S0_filename == NULL) {
            fprintf(stderr, "Option 'S0compute' requires output filename to be\n");
            fprintf(stderr, "specified using -S0 is you want to save it.\n");
        } else {
            fprintf(stderr, "S0 image will be computed using internal methods and\n");
            fprintf(stderr, "will be saved to: %s.\n", mow->S0_filename);
        }
    } else if (mow->S0compute == 0) {
        if (mow->S0_filename == NULL) {
            fprintf(stderr, "S0 image not specified. Will compute using internal method.\n");
            fprintf(stderr, "Resulting S0 computation will not be saved.\n");
            mow->S0compute = 1;
        } else {
            fprintf(stderr, "S0 image will be loaded from %s.\n", mow->S0_filename);
            mow->diff->S0 = nifti_image_read(mow->S0_filename, 1);
            if (mow->diff->S0 == NULL) {
                exit(1);
            }
        }
    }
    
    if (mow->data_filename == NULL) {
        fprintf(stderr, "The -data <path> option is required to specify the location\n");
        fprintf(stderr, "of the diffusion-weighted data.\n");
        exit(1);
    } else if (mow->bval_filename == NULL) {
        fprintf(stderr, "The -bval <path> option is required to specify the location\n");
        fprintf(stderr, "of the b-values for this data set.\n");
        exit(1);
    } else if (mow->bvec_filename == NULL) {
        fprintf(stderr, "The -bvec <path> option is required to specify the location\n");
        fprintf(stderr, "of the gradient directions for this data set.\n");
        exit(1);
    }
    
    /* these functions call exit() if anything bad happens. */
    fprintf(stderr, "Loading diffusion data...\n");
    read_diff_data_from_file(mow->data_filename, mow->diff);
    fprintf(stderr, "Loading bvalues...\n");
    read_bvals_from_file(    mow->bval_filename, mow->diff);
    fprintf(stderr, "Loading gradient vectors...\n");
    read_bvecs_from_file(    mow->bvec_filename, mow->diff);
    
    fprintf(stderr, "Loading binary mask...\n");
    mow->diff->mask = nifti_image_read(mow->mask_filename, 1);
    if (mow->diff->mask == NULL) {
        fprintf(stderr, "Unable to load mask file: %s\n", mow->mask_filename);
        fprintf(stderr, "I will attempt to reconstruct ALL voxels, even those\n");
        fprintf(stderr, "that may lie outside of the subject matter. Use results\n");
        fprintf(stderr, "with caution.\n\n");
        //exit(1);
    } else {
        nii_recast_to_int32(mow->diff->mask);    
    }
    
}

/*****************************************************************************
 
 *****************************************************************************/
void mow_print_usage(void)
{
    
    fprintf(stderr, "Usage: mow_recon [OPTIONS]\n\n");
    fprintf(stderr, "where OPTIONS is all of the following:\n");
    fprintf(stderr, "  -data <path>      path to diffusion-weighted data in\n");
    fprintf(stderr, "                    NIFTI format (.nii or .nii.gz)\n\n");
    fprintf(stderr, "  -bval <path>      path to the text file containing the\n");
    fprintf(stderr, "                    b-values corresponding to the diffusion\n");
    fprintf(stderr, "                    weighted volumes in the data file.\n\n");
    fprintf(stderr, "  -bvec <path>      path to the text file containing the\n");
    fprintf(stderr, "                    gradient directions for the diffusion\n");
    fprintf(stderr, "                    weighted volumes in the data file.\n\n");
    fprintf(stderr, "  -mask <path>      path to a binary mask with nonzero entries\n");
    fprintf(stderr, "                    indicating which voxels should be considered.\n\n");
    fprintf(stderr, "  -S0   <path>      path to an image containing the S0 data.\n");
    fprintf(stderr, "                    This is optional, and if not supplied the\n");
    fprintf(stderr, "                    S0 data will be computed internallly.\n\n");
    fprintf(stderr, "  -ndir <integer>   Specifies the number of directions to save\n");
    fprintf(stderr, "                    per voxel. Maxima will be sorted by probability\n");
    fprintf(stderr, "                    and saved from highest probability to lowest.\n\n");
    fprintf(stderr, "  -pthresh <float>  Specifies the minimum probability (from 0 to 1)\n");
    fprintf(stderr, "                    that a maximal direction must have in order to be saved.\n\n");
    fprintf(stderr, "  -delta_lg <float> Specifies the 'large delta' which indicates the time\n");
    fprintf(stderr, "                    between diffusion pulsed in the pulse sequence.\n\n");
    fprintf(stderr, "  -delta_sm <float> Specifies the 'small delta' which indicates the time that\n");
    fprintf(stderr, "                    the diffusion pulse is turned on.\n\n");
    fprintf(stderr, "  -radius <float>   The displacement radius to consider.\n\n");
    fprintf(stderr, "  -odir <dirpath>   Path to a directory where the output files\n");
    fprintf(stderr, "                    should be saved. They will be named 'V_xx_all.nii\n");
    fprintf(stderr, "                    where xx ranges from 0 to ndir-1.\n\n");
    fprintf(stderr, "  -datadir <path>   path to the TrackTools data directory.\n\n");
    fprintf(stderr, "  -log-bad-voxels   If this option is used, mow_recon will print\n");
    fprintf(stderr, "                    a list of all voxels that could not be reconstructed.\n\n");
    exit(0);
    
}

/*****************************************************************************
 *  Initializes the OUTPUT_DATA struct based on the dimensions and orientation
 *  of the template. It returns an properly allocated OUTPUT_DATA struct
 *  which can be passed to add_maxima_to_output.
 *
 *  It allocates space for nfiles.
 *****************************************************************************/
OUTPUT_DATA *initialize_output (nifti_image *template, int nfiles)
{
    int i;
    
    nifti_image      **out = malloc(sizeof(nifti_image *)      * nfiles);
    nifti_brick_list **nbl = malloc(sizeof(nifti_brick_list *) * nfiles);
    if (out == NULL || nbl == NULL) {
        fprintf(stderr, "Unable to allocate memory for nifti story in output data.\n");
        return NULL;
    }
    
    for (i=0 ; i<nfiles; i++) {
        
        out[i] = nifti_simple_init_nim();
        if (out[i] == NULL) {
            fprintf(stderr, "Unable to allocate memory for nifti image.\n");
            return NULL;
        }
        
        memcpy(out[i], template, sizeof(*out[i]));
        
        out[i]->datatype  = DT_FLOAT32;
        out[i]->nbyper    = 4;
        out[i]->nt        = 3;
        out[i]->nvox      = out[i]->nx*out[i]->ny*out[i]->nz*out[i]->nt;
        out[i]->dim[4]    = 3;
        out[i]->fname     = NULL;
        out[i]->iname     = NULL;
        out[i]->scl_slope = 1.0;
        out[i]->scl_inter = 0.0;
        out[i]->cal_max   = 0.0;
        out[i]->cal_min   = 0.0;
        out[i]->nifti_type= NIFTI_FTYPE_NIFTI1_1;	
        sprintf(out[i]->descrip, "TrackTools MOW Reconstruction");
        
        nbl[i] = malloc(sizeof(nifti_brick_list));
        if (nbl[i] == NULL) {
            fprintf(stderr, "Unable to allocate memory for nifti brick list in output data.\n");
            return NULL;
        }
        
        nbl[i]->nbricks = 3;
        nbl[i]->bsize   = out[i]->nx*out[i]->ny*out[i]->nz*out[i]->nbyper;
        nbl[i]->bricks  = malloc(sizeof(void *) * 3);
        if (nbl[i]->bricks == NULL) {
            fprintf(stderr, "Unable to allocate memory for data storage in output data.\n");
            return NULL;
        }
        
        nbl[i]->bricks[0] = malloc(nbl[i]->bsize);
        nbl[i]->bricks[1] = malloc(nbl[i]->bsize);
        nbl[i]->bricks[2] = malloc(nbl[i]->bsize);
        if (nbl[i]->bricks[0] == NULL ||
            nbl[i]->bricks[1] == NULL ||
            nbl[i]->bricks[2] == NULL) {
            fprintf(stderr, "Unable to allocate nifti bricks for output data.\n");
            return NULL;
        }
    }
    
    OUTPUT_DATA *output = malloc(sizeof(OUTPUT_DATA));
    if (output == NULL) {
        fprintf(stderr, "Unable to allocate memory for output data.\n");
        return NULL;
    }
    
    output->nim = out;
    output->nbl = nbl;
    output->num_images = nfiles;
    
    return output;
    
}

/*****************************************************************************
 * Adds the directions specified in maxima_list to the output data structure
 * at voxel location (x,y,z). Note, each maxima coordinate is added to a brick
 * in the brick list. That is, (x -> brick 0, y -> brick 1, z-> brick 2) and
 * the process is repeated for each maxima.
 *****************************************************************************/
void add_maxima_to_output(OUTPUT_DATA *output, int x, int y, int z, 
                          float **vertlist, MAXIMA *maxima_list, int n_maxima)
{
    
    int i;
    int index = nii_voxel3_index(output->nim[0], x, y, z);
    
    if (n_maxima > output->num_images) {
        /* fprintf(stderr, "  WARNING: Voxel: [%d,%d,%d], adding more maxima than number of files (%d > %d).\n",
         x, y, z, n_maxima, output->num_images); */
        n_maxima = output->num_images;
    }
    
    for (i=0; i<n_maxima; i++) {
        float *mvert = vertlist[ maxima_list[i].index ];
        ((float *)output->nbl[i]->bricks[0])[index] = mvert[0]; 
        ((float *)output->nbl[i]->bricks[1])[index] = mvert[1]; 
        ((float *)output->nbl[i]->bricks[2])[index] = mvert[2]; 
    }
    
}

/*****************************************************************************
 * Writes the OUTPUT_DATA structure to disk. One file for each nifti_image
 * in the OUTPUT_DATA structure. 
 *
 * Files are saved in basedir, and are strictly named to be compatible with
 * track_tracker.
 *
 * **NOTE** existing files will be overwritten quietly.
 *****************************************************************************/
int save_output (const char *basedir, OUTPUT_DATA *output)
{
    int i;
    char *output_iname = malloc(sizeof(char)*1024);
    char *output_fname = malloc(sizeof(char)*1024);
    if (output_iname == NULL || output_fname == NULL) {
        fprintf(stderr, "Unable to save output. Out of memory.\n");
        exit(1);
    }
    
    if (strlen(basedir)+13 > (1024-1)) {
        fprintf(stderr, "Output file name too long.\n");
        return 0;
    }
    
    for (i=0; i<output->num_images; i++) {
        sprintf(output_iname, "%s%cV_%02d_all.nii", basedir, DIRSEP, i);
        sprintf(output_fname, "%s%cV_%02d_all.nii", basedir, DIRSEP, i);
        znzFile fp = znzopen(output_fname, "wb", 0);
        if (fp != NULL) {
            output->nim[i]->iname = output_iname;
            output->nim[i]->fname = output_fname;
            nifti_image_write_hdr_img2(output->nim[i], 1, "wb", fp, output->nbl[i]);
        } else {
            fprintf(stderr, "Unable to open output file for writing (%s).\n",
                    output_fname);
            exit(1);
        }
    }
    
    return 1;
    
}

/*****************************************************************************
 * Reads a text file containing the bvalues into a floating point array
 * contained in the DIFF_DATA structure.
 *
 * Also separates any low-bvalues from high-bvalues according to their indices
 * in the floating point array.
 *****************************************************************************/
void read_bvals_from_file(char *filename, DIFF_DATA *diff)
{
    
    int i, n_low, n_high;
    float *bval_array = malloc(sizeof(float) * diff->n_volumes);
    if (bval_array == NULL) {
        fprintf(stderr, "Unable to allocate memory for bvalue table.\n");
        exit(1);
    }
    
    int n_bvals = read_acsii_file_to_float_array(filename, bval_array,
                                                 diff->n_volumes);
    diff->bvals = bval_array;
    
    n_low = 0;
    for (i=0; i<n_bvals; i++) {
        if (bval_array[i] < 200.0) {
            n_low++;
        }
    }
    
    diff->b_low_ind = malloc(sizeof(int) * n_low);
    if (diff->b_low_ind == NULL) {
        fprintf(stderr, "Unable to allocate memory for low-b data.\n");
        exit(1);
    }
    
    diff->b_high_ind = malloc(sizeof(int) * n_bvals-n_low);
    if (diff->b_high_ind == NULL) {
        fprintf(stderr, "Unable to allocate memory for high-b data.\n");
        exit(1);
    }
    
    n_low = 0;
    n_high = 0;
    for (i=0; i<n_bvals; i++) {
        if (bval_array[i] < 200.0) {
            diff->b_low_ind[n_low++] = i;
        } else {
            diff->b_high_ind[n_high++] = i;
        }
    }
    
    diff->n_b_low = n_low;
    diff->n_b_high = n_high;
    
}

/*****************************************************************************
 * Reads a text file containing gradient directions.
 * Text file should be organized this way:
 *
 *   for n = 1 ... number of diffusion scans:
 *   first line: g1_x ... gn_x
 *   next line : g1_y ... gn_y
 *  third line : g1_z ... gn_z
 *
 *****************************************************************************/
void read_bvecs_from_file(char *filename, DIFF_DATA *diff)
{
    
    float *bvec_array = malloc(sizeof(float) * 3*diff->n_volumes);
    if (bvec_array == NULL) {
        fprintf(stderr, "Unable to allocate memory for gradient table.\n");
        exit(1);
    }
    
    int n_bvecs = read_acsii_file_to_float_array(filename, bvec_array,
                                                 3*diff->n_volumes);
    
    if (n_bvecs != 3*diff->n_volumes) {
        fprintf(stderr, "Something's wrong with the bvecs: expected: %d read: %d\n",
                3*diff->n_volumes, n_bvecs);
        exit(1);
    }
    
    diff->bvecs = bvec_array;
    
}

/*****************************************************************************
 * Loads diffusion-weighted data from a nifti file.
 *****************************************************************************/
void read_diff_data_from_file(char *filename, DIFF_DATA *diff)
{
    
    nifti_brick_list *nbl = NULL;
    nifti_image *nim = NULL;
    
    nbl = malloc(sizeof(nifti_brick_list));
    if (nbl == NULL) {
        fprintf(stderr, "Unable to allocate memory for nifti data.\n");
        exit(1);
    }
    
    nim = nifti_image_read_bricks(filename, 0, NULL, nbl);
    if (nim == NULL) { 
        fprintf(stderr, "Unable to read diffusion data file.\n");
        exit(1);
    }
    
    /* we use scl_slope, and if == 0, then should be disregarded. */
    if (abs(nim->scl_slope) < 1e-5) {
        nim->scl_slope = 1.0;
    }
    
    diff->nii_brick_list = nbl;
    diff->nii_image      = nim;
    diff->n_volumes      = nbl->nbricks;
    
    diff->single_voxel_storage = malloc(sizeof(double) * nbl->nbricks);    
    if (diff->single_voxel_storage == NULL) {
        fprintf(stderr, "Unable to allocate memory for diffusion voxel storage.\n");
        exit(1);
    }
}

/*****************************************************************************
 * Reads a text file into a floating point array. Note that file is assumed to 
 * be a text file and behavior is undefined if it is not. 
 *****************************************************************************/
int read_acsii_file_to_float_array(char *filename, float *data, int data_size)
{
    
    /* I may have to make this smaller depending on the stack size */
    char *strbuf = malloc(sizeof(char) * 10240+1);
    float *tmp   = malloc(sizeof(float) * 10240+1);
    
    int elements_read = 0;
    
    FILE *f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "Unable to open bval file: %s\n", filename);
        perror("Reason");
        return 0;
    }
    
    while (! feof(f)) {
        int len = 0;
        int num_toks = 0;
        int i=0;
        
        memset(strbuf, 0, sizeof(char)*10240);
        
        fgets(strbuf, 10239, f);
        len = strlen(strbuf);
        if (len <= 1) continue;
        
        strbuf[len-1] = '\0';
        
        __tokenize_line(strbuf, " \t", tmp, &num_toks);
        
        for (i=0; i<num_toks && elements_read+i < data_size; i++) {
            data[elements_read+i] = tmp[i];
        }
        elements_read += i;
        
        if (elements_read >= data_size) {
            break;
        }
        
    }
    
    fclose(f);
    free(strbuf);
    free(tmp);
    
    return elements_read;
}

/*****************************************************************************
 * Uses the find the relative maxima of the spherical function over the 
 * spherical domain. Note that each domain file contains a connectivity list
 * that allows one to find the immediate neighbors of a vertex in O(1) time.
 *
 * Now, *restart contains a spherically distributed subset of *domain, so that
 * by started at each vertex in *restart, each local maxima can be found.
 *
 * Finally, duplicates and antipodally symmetric vertices are removed, and the
 * maxima are sorted by value.
 * 
 * returns the number of vertices found.
 *****************************************************************************/

/* size of *values should be equal to number of elements in *domain->vertices */
int find_local_maxima(ICOS_TESS *domain, double *values, double min_value,
                      ICOS_TESS *restart, MAXIMA *maxima_list) 
{
    
    int v1, v2, j;
    int n_maxima = 0;
    float dp;
    /* MAXIMA *maxima_list; */
    
    /* these are the indices into the array of vertices. */
    for (v1=0; v1<restart->num_vertices; v1++) {
        
        int rs_vert = v1; //restart->vert[v1];
        
        double curval = values[rs_vert]; /* the value of the function a
                                          vertex[rs_vert] */
        int stop = 0;
        
        while (stop == 0) {
            
            int maxind = rs_vert;
            double maxval = 0;
            
            /* the number of neighbors of this vertex */
            int n_neighbors = domain->connectivity[rs_vert][0];
            
            /* find the neighbor that has the largest value greater than
             the current value */
            for(v2=1; v2<=n_neighbors; v2++) {
                
                /* value of the function at this particular neighbor */
                double testval = values[domain->connectivity[rs_vert][v2]];
                
                if (testval > maxval) {
                    maxval = testval;
                    maxind = domain->connectivity[rs_vert][v2];
                }
                
            }
            
            if (maxval > curval) {
                /* found better vertex */
                curval  = maxval;
                rs_vert = maxind;
            } else {
                /* no neighbor has higher value; we are at a maxima */
                stop = 1;
            }
            
        }
        
        /* rs_vert is the index into the vertex array of the vertex that holds
         the maxima.
         curval is the value of the funcation a vertices[rs_vert].
         */
        
        for (j=0; j<n_maxima; j++) {
            float *pb_vert = domain->vertices[rs_vert];
            float *cp_vert = domain->vertices[maxima_list[j].index];
            dp = fabs(pb_vert[0]*cp_vert[0] +
                      pb_vert[1]*cp_vert[1] + 
                      pb_vert[2]*cp_vert[2] );
            if (dp > 0.939693) {
                /* reject this maxima, duplicate */
                break;
            }
        }
        
        if (j == n_maxima && values[rs_vert] >= min_value) {
            /* keep this maxima */
            maxima_list[n_maxima].index = rs_vert;
            maxima_list[n_maxima].value = values[rs_vert];
            n_maxima++;
        }
        
        /*
         printf("[ %f ][ %f, %f, %f ]\n", curval,
         domain->vertices[rs_vert][0],
         domain->vertices[rs_vert][1],
         domain->vertices[rs_vert][2]);
         */
        
    }
    
    qsort((void *)maxima_list, n_maxima, sizeof(MAXIMA), &maxima_compare);
    
    return n_maxima;
    
}

/*****************************************************************************
 * This is not implemented. NNLS is used in the main loop to compute this.
 *****************************************************************************/

double *compute_reconstruction_weights(MOW_RECON *mow, int x, int y, int z)
{
    double *d = NULL;
    return d;
}

/*****************************************************************************
 * Builds an array of B matrices using the bvals and gradient directions.
 * Only the upper triangular entries are stored. Note that the off-diagonal
 * entries are multiplied by 2. This is a MAS convention that makes computing
 * DTI easier. We do it here, and note that the factor of two is removed when
 * the A matrix is created. 
 *****************************************************************************/

void compute_bmatrix(DIFF_DATA *diff)
{
    int i;
    
    float *bvec = diff->bvecs;
    
    diff->b_matrices = malloc(sizeof(SYMMAT33 *) * diff->n_volumes);
    if (diff->b_matrices == NULL) {
        fprintf(stderr, "Unable to allocate memory for b matrices.\n");
        exit(1);
    }
    
    for (i=0; i<diff->n_volumes; i++) {
        float gradvec[3];
        gradvec[0] = bvec[diff->n_volumes*0+i];
        gradvec[1] = bvec[diff->n_volumes*1+i];
        gradvec[2] = bvec[diff->n_volumes*2+i];
        float bval = diff->bvals[i]*1.0e-3;
        SYMMAT33 *bmatrix = malloc(sizeof(SYMMAT33));
        if (bmatrix == NULL) {
            fprintf(stderr, "Unable to allocate b-matrix.\n");
            exit(1);
        }
        
        bmatrix->data[0] =     bval * gradvec[0] * gradvec[0]; /* xx */
        bmatrix->data[1] =     bval * gradvec[1] * gradvec[1]; /* yy */
        bmatrix->data[2] =     bval * gradvec[2] * gradvec[2]; /* zz */
        bmatrix->data[3] = 2 * bval * gradvec[0] * gradvec[1]; /* xy */
        bmatrix->data[4] = 2 * bval * gradvec[0] * gradvec[2]; /* xz */
        bmatrix->data[5] = 2 * bval * gradvec[1] * gradvec[2]; /* yz */
        
        diff->b_matrices[i] = bmatrix;
    }
    
}

/*****************************************************************************
 * This reconstruction matrix is computed using Eqn 16 from the paper.
 *****************************************************************************/

double *make_recon_matrix(MOW_RECON *mow)
{
    
    ICOS_TESS *deco_tess = mow->deco_tess;
    ICOS_TESS *reco_tess = mow->reco_tess;
    
    double *reco_matrix;
    double diff_rad_sq = mow->diff_radius * mow->diff_radius;
    double diff_time   = mow->diff->delta_lg - mow->diff->delta_sm/3.0; 
    
    int n_deco_dirs = deco_tess->num_vertices;
    int n_reco_dirs = reco_tess->num_vertices;
    
    int i,j;
    
    reco_matrix = malloc(sizeof(double)*n_reco_dirs*n_deco_dirs);
    if (reco_matrix == NULL) {
        fprintf(stderr, "Unable to allocate memory for reconstruction matrix\n");
        return NULL;
    }
    
    for (i=0; i<n_deco_dirs; i++) {
        
        for (j=0; j<n_reco_dirs; j++) {
            
            float *r = reco_tess->vertices[j];
            double rt[6];
            double rt_d_r = 0;
            double exp_arg;
            int n;
            
            rt[0] = r[0]*r[0] * diff_rad_sq;
            rt[1] = r[1]*r[1] * diff_rad_sq;
            rt[2] = r[2]*r[2] * diff_rad_sq;
            rt[3] = r[2]*r[1] * diff_rad_sq;
            rt[4] = r[0]*r[1] * diff_rad_sq;
            rt[5] = r[0]*r[2] * diff_rad_sq;
            
            for (n=0; n<6; n++) 
                rt_d_r += mow->D_i_inv[i]->data[n] * (double)rt[n];
            
            exp_arg = -(rt_d_r/(4.0 * diff_time));
            
            reco_matrix[n_reco_dirs*i + j] = exp(exp_arg);
            
        }
        
    }
    
    mow->reco_matrix = reco_matrix;
    
    return reco_matrix;
    
}

/*****************************************************************************
 * Creates the A matrix from the paper. Page 168, including Eqn #15
 *****************************************************************************/

double *make_A_matrix(MOW_RECON *mow) 
{
    
    MAT33 *tx1 = NULL;
    MAT33 *tx2 = NULL;
    MAT33 *tx = NULL;
    MAT33 *Q_mat = NULL, *Q_mat_t = NULL, *diag_evals = NULL;
    MAT33 *Sig_i = NULL, *D_i_inverse = NULL, *mat_scratch = NULL;
    
    int i,j;
    float pi_frac = 3.0;
    float arg = M_PI/pi_frac;
    float sa = sin(arg);
    float ca = cos(arg);
    
    ICOS_TESS *tess = mow->deco_tess;
    DIFF_DATA *diff = mow->diff;
    double *A_matrix;
    
    SYMMAT33 **b_matrix = diff->b_matrices;
    
    /* transformation matrix used to obtain orthogonal vector set */
    tx = malloc(sizeof(MAT33));
    
    /* required for mow calculation */
    Q_mat       = malloc(sizeof(MAT33));
    Q_mat_t     = malloc(sizeof(MAT33));
    diag_evals  = malloc(sizeof(MAT33));
    
    /* this is temporary placeholder used for multiplying three matrices
     together */
    mat_scratch = malloc(sizeof(MAT33));
    
    if (tx == NULL || Q_mat == NULL || Q_mat_t == NULL || diag_evals == NULL ||
        mat_scratch == NULL) {
        fprintf(stderr, "Unable to allocate memory for computation matrices.\n");
        return NULL;
    }
    
    A_matrix = malloc(sizeof(double) * tess->num_vertices * diff->n_b_high);
    if (A_matrix == NULL) {
        fprintf(stderr, "Unable to allocate memory for A matrix.\n");
        return NULL;
    }
    
    mow->D_i_inv = malloc(sizeof(SYMMAT33 *) * tess->num_vertices);
    if (mow->D_i_inv == NULL) {
        fprintf(stderr, "Unable to allocate memory for D_inverse.\n");
        return NULL;
    }
    
    /* matrices used to create an orthogonal set of vectors with the deconvolution
     vector as primary eigenvector. The secondary and tertiary evectors can be
     oriented any way. These matrices rotate the original deconvolution vector 
     out of the plane, then the second evector is created by the cross product
     of the rotated vector with the original. Then, the final evector is created
     by the cross product of the original vector and the secondary. */
    tx1 = MAT33_make(ca, -sa , 0.0, 
                     sa,  ca , 0.0,
                     0.0, 0.0, 1.0);
    
    tx2 = MAT33_make(ca , 0.0, sa, 
                     0.0, 1.0, 0.0,
                     -sa, 0.0, ca);
    
    MAT33_mult(tx2, tx1, tx);
    
    /* I am just going to reuse the memory of tx1 and tx2 since we don't need
     them anymore */
    Sig_i       = tx1;
    D_i_inverse = tx2;
    
    for (i=0; i<tess->num_vertices; i++) {
        float tmp[3];
        float egv1[3], egv2[3];
        float *tmp_v = tess->vertices[i];
        
        /* this is tx * tmp_v matrix mult. */
        tmp[0] = tmp_v[0]*tx->data[0] + tmp_v[1]*tx->data[1] + tmp_v[2]*tx->data[2];
        tmp[1] = tmp_v[0]*tx->data[3] + tmp_v[1]*tx->data[4] + tmp_v[2]*tx->data[5];
        tmp[2] = tmp_v[0]*tx->data[6] + tmp_v[1]*tx->data[7] + tmp_v[2]*tx->data[8];
        VEC_make_unit(tmp);
        
        /* tmp is a preturbed version of the original */
        VEC_crossp(tmp, tmp_v, egv1);
        VEC_make_unit(egv1);
        
        /* evec 2 created from cross product */
        VEC_crossp(egv1, tmp_v, egv2);
        VEC_make_unit(egv2);
        
        MAT33_assign(tmp_v[0], egv1[0], egv2[0], 
                     tmp_v[1], egv1[1], egv2[1],
                     tmp_v[2], egv1[2], egv2[2], Q_mat);
        
        MAT33_assign(tmp_v[0], tmp_v[1], tmp_v[2],
                     egv1[0],  egv1[1],  egv1[2],
                     egv2[0],  egv2[1],  egv2[2], Q_mat_t);
        
        MAT33_assign(mow->deco_evals[0], 0, 0,
                     0, mow->deco_evals[1], 0,
                     0, 0, mow->deco_evals[2], diag_evals);
        
        /* Q_mat * diag_evals * Q_mat_transpose */
        /* MAT33_mult( MAT33_mult(Q_mat_t, diag_evals, mat_scratch) , Q_mat, Sig_i); */
        MAT33_mult(Q_mat_t, diag_evals, mat_scratch);
        MAT33_mult(mat_scratch, Q_mat, Sig_i);
        
        /* now we invert diag_evals easily since it is diag. */
        for (j=0; j<9; j++) {
            if (diag_evals->data[j] != 0.0)
                diag_evals->data[j] = 1.0/diag_evals->data[j];
        }
        
        /* D_I = Sig_i before dividing by p, so D_i_inverse should
         happen here (diag_evals has already been inverted). */
        MAT33_mult(Q_mat_t, diag_evals, mat_scratch);
        MAT33_mult(mat_scratch, Q_mat, D_i_inverse);
        
        /* D_i_inv gets saved for reconstrucion later */
        mow->D_i_inv[i] = malloc(sizeof(SYMMAT33));
        /* go ahead and compute elements for quadratic form */
        mow->D_i_inv[i]->data[0] = D_i_inverse->data[0];
        mow->D_i_inv[i]->data[1] = D_i_inverse->data[4];
        mow->D_i_inv[i]->data[2] = D_i_inverse->data[8];
        mow->D_i_inv[i]->data[3] = 2*D_i_inverse->data[5];
        mow->D_i_inv[i]->data[4] = 2*D_i_inverse->data[1];
        mow->D_i_inv[i]->data[5] = 2*D_i_inverse->data[2];
        /* create Sig_i/p */
        for (j=0; j<9; j++) Sig_i->data[j] /= DECO_P;
        
        /* build A as per eqn XXX in the paper */
        for (j=0; j<diff->n_b_high; j++) {
            
            int ind = diff->b_high_ind[j];
            float trace = 0;
            /* we take out the factor of two on the off-diagonals here */
            mat_scratch->data[0] = b_matrix[ind]->data[0];
            mat_scratch->data[1] = b_matrix[ind]->data[3]/2.0;
            mat_scratch->data[2] = b_matrix[ind]->data[4]/2.0;
            
            mat_scratch->data[3] = b_matrix[ind]->data[3]/2.0;
            mat_scratch->data[4] = b_matrix[ind]->data[1];
            mat_scratch->data[5] = b_matrix[ind]->data[5]/2.0;
            
            mat_scratch->data[6] = b_matrix[ind]->data[4]/2.0;
            mat_scratch->data[7] = b_matrix[ind]->data[5]/2.0;
            mat_scratch->data[8] = b_matrix[ind]->data[2];
            
            MAT33_mult(Sig_i, mat_scratch, Q_mat);
            trace = Q_mat->data[0] + Q_mat->data[4] + Q_mat->data[8];
            A_matrix[diff->n_b_high*i + j] = (double) pow(1.0 + trace, -DECO_P);
            
        }
    }
    
    replace_realsmall_w_zeros(A_matrix, tess->num_vertices * diff->n_b_high,
                              (double) 1.0e-10);
    
    mow->A_matrix = A_matrix;
    
    free(tx);
    free(Sig_i);
    free(D_i_inverse);
    free(Q_mat);
    free(Q_mat_t);
    free(diag_evals);
    free(mat_scratch);
    
    return A_matrix;
}

/*****************************************************************************
 * Loads a spherical domain into the appropiate data structure. See 
 * mow_recon.h for a description of the data structure. Returns NULL if
 * memory allocation fails or some other IO problem prevents the data from
 * loading.
 *****************************************************************************/

ICOS_TESS *load_tess_from_file(const char *tess_file)
{
    
    int i, vert, iobytes;
    int num_dimensions;
    ICOS_TESS *tess;
    
    FILE *f = fopen(tess_file, "rb");
    if (f == NULL) {
        fprintf(stderr, "Unable to open tessellation file: %s\n", tess_file);
        perror("Reason");
        return NULL;
    }
    
    tess = malloc(sizeof(ICOS_TESS));
    if (tess == NULL) {
        fprintf(stderr, "Unable to allocate memory for tessellation.\n");
        return NULL;
    }
    
    iobytes = fread(&num_dimensions, 1, 4, f);
    iobytes += fread(&(tess->num_vertices), 1, 4, f);
    if (iobytes != 8) {
        fprintf(stderr, "Error reading from tessellation file (1).\n");
        perror("Reason");
        fclose(f);
        return NULL;
    }
    
    tess->vertices = malloc(sizeof(float *) * tess->num_vertices);
    if (tess->vertices == NULL) {
        fprintf(stderr, "Unable to allocate memory for vertex array.\n");
        return NULL;
    }
    
    for (i=0; i<tess->num_vertices; i++) {
        tess->vertices[i] = malloc(sizeof(float) * 3);
        if (tess->vertices[i] == NULL) {
            fprintf(stderr, "Unable to allocate memory for tess. vertex.\n");
            return NULL;
        }
        
        iobytes = fread(tess->vertices[i], 3, 4, f);
        if (iobytes != 4) {
            fprintf(stderr, "Unable to read from tessellation file (2) (%d).\n",
                    iobytes);
            perror("Reason");
            fclose(f);
            return NULL;
        }
    }
    
    tess->connectivity = malloc(sizeof(int *) * tess->num_vertices);
    if (tess->connectivity == NULL) {
        fprintf(stderr, "Unable to allocate memory for connectivity array.\n");
        return NULL;
    }
    
    for (vert=0; vert<tess->num_vertices; vert++) {
        int conn_size = 0;
        
        iobytes = fread(&conn_size, 1, 4, f);
        if (iobytes != 4) {
            fprintf(stderr, "Error reading connectivity from file.\n");
            perror("Reason");
            fclose(f);
            return NULL;
        }
        tess->connectivity[vert] = malloc(sizeof(int)*(conn_size+1));
        if (tess->connectivity[vert] == NULL) {
            fprintf(stderr, "Unable to allocate memory for connectivity list\n");
            return NULL;
        }
        
        tess->connectivity[vert][0] = conn_size;
        
        iobytes = fread(tess->connectivity[vert]+1, conn_size, 4, f);	
        if (iobytes != 4) {
            fprintf(stderr, "Error reading connectivity from file.\n");
            perror("Reason");
            fclose(f);
            return NULL;
        }
        
    }
    
    fclose(f);
    
    return tess;
    
}

/*****************************************************************************
 * Internal function to turn a string into an array of floats. Used in reading
 * bval/bvec files.
 *****************************************************************************/

void __tokenize_line(char *line, const char *tokens, float *split, int *nsplit)
{
    
    int n = 0;
    char *tmpstr = strtok(line, " \t");
    split[n++] = (float) atof(tmpstr);
    
    while ( NULL != (tmpstr = strtok(NULL, " \t\n")) ) {
        
        if (strcmp(tmpstr, "\n") == 0) continue;
        
        split[n++] = (float) atof(tmpstr);
        
        /* printf("%d, %d, %+.2f\n", i, (int)strlen(tmpstr), (float)atof(tmpstr)); */
    }
    
    *nsplit = n;
}

/*****************************************************************************
 * Function used by qsort to compare two maxima.
 * Comparison is based on probability values.
 *****************************************************************************/

int maxima_compare(const void *m1, const void *m2)
{
    
    MAXIMA *__m1 = (MAXIMA *) m1;
    MAXIMA *__m2 = (MAXIMA *) m2;
    
    if (__m1->value == __m2->value) {
        return 0;
    } else {
        return (__m1->value < __m2->value) ? 1 : -1;
    }
    
}

/*****************************************************************************
 * Function to replace small double values with zero.
 *****************************************************************************/

void replace_realsmall_w_zeros(double *input, int size, double tol)
{
    int i;
    
    for (i=0; i<size; i++) {
        if (fabs(input[i]) < tol)
            input[i] = 0.0;
    }
    
}

/*****************************************************************************
 * Loads a voxel's diffusion data and divides it by the voxel's S0. 
 * - Tests to see if the voxel is masked by diff->mask.
 * - Returns only the diffusion-weighted values (b > 200)
 * - returns 1 is loading succeeds
 * - returns -1 if a inf/nan situation prevents the data (ie S0 = 0)
 * - returns 0 is a mask hit prevents the data from loading.
 *****************************************************************************/

int load_voxel_double_highb(DIFF_DATA *diff, int x, int y, int z) {
    
    int index = nii_voxel3_index(diff->nii_image, x, y, z);
    int i;
    double *data = diff->single_voxel_storage;
    double scl_slope = (double) diff->nii_image->scl_slope;
    double scl_inter = (double) diff->nii_image->scl_inter;
    float S0;
    
    if (diff->mask != NULL) {
        int *mask = (int *)diff->mask->data;
        if (mask != NULL && mask[index] == 0) return 0;
    }
    
    S0 = (double) read_nii_voxel_anytype(diff->S0->data, 
					 index, 
					 diff->S0->datatype);
    if (S0 == 0) return -1;
    
    nifti_brick_list *nbl = diff->nii_brick_list;
    
    for (i=0; i<diff->n_b_high; i++) {
        
        int b_ind = diff->b_high_ind[i];
	
	data[i] = (double)read_nii_voxel_anytype(nbl->bricks[b_ind], 
						 index, 
						 diff->nii_image->datatype);
	
        if (isinf(data[i]) || isnan(data[i])) {
            return -2;
        }
        
        data[i] = (data[i]*scl_slope+scl_inter)/S0; 
        
    }
    
    return 1;
    
}

/*****************************************************************************
 * Loads a voxel's diffusion data and DOES NOT divide it by the voxel's S0. 
 * - Tests to see if the voxel is masked by diff->mask.
 * - Returns all data.
 * - returns 1 is loading succeeds
 * - returns -1 if a inf/nan situation prevents the data (ie S0 = 0)
 * - returns 0 is a mask hit prevents the data from loading.
 *****************************************************************************/

int load_voxel_double_all(DIFF_DATA *diff, int x, int y, int z, double *dest) {
    
    int index = nii_voxel3_index(diff->nii_image, x, y, z);
    int i;
    double *data;
    double scl_slope = (double) diff->nii_image->scl_slope;
    double scl_inter = (double) diff->nii_image->scl_inter;
    
    if (dest == NULL) {
        data = diff->single_voxel_storage;
    } else {
        data = dest;
    }
    
    if (diff->mask != NULL) {
        int *mask = (int *)diff->mask->data;
        if (mask != NULL && mask[index] == 0) return 0;
    }
    
    nifti_brick_list *nbl = diff->nii_brick_list;
    
    for (i=0; i<diff->n_volumes; i++) {

	data[i] = (double)read_nii_voxel_anytype(nbl->bricks[i], 
						 index, 
						 diff->nii_image->datatype);

        if (isinf(data[i]) || isnan(data[i])) {
            return -1;
        }
        
        data[i] = data[i]*scl_slope+scl_inter;
    }
    
    return 1;
    
}

/*****************************************************************************
 * Loads a voxel's data and casts from any type to double. 
 * returns 0 on undefined type. 
 *****************************************************************************/

double read_nii_voxel_anytype(void *src, int index, int datatype)
{
    double dest = 0;
    
    switch (datatype) {
	case DT_UINT8:
	    dest = (double) ( (unsigned char *)(src))[index];
	    break;
	case DT_INT8:
	    dest = (double) (          (char *)(src))[index];
	    break;
	case DT_UINT16:
	    dest = (double) ((unsigned short *)(src))[index];
	    break;
	case DT_INT16:
	    dest = (double) (         (short *)(src))[index];
	    break;
	case DT_UINT32:
	    dest = (double) (  (unsigned int *)(src))[index];
	    break;
	case DT_INT32:
	    dest = (double) (           (int *)(src))[index];
	    break;
	case DT_FLOAT32:
	    dest = (double) (         (float *)(src))[index];
	    break;
	case DT_FLOAT64:
	    dest = (                 (double *)(src))[index];
	    break;
	default: 
	    fprintf(stderr, "Nifti datatype not supported: %s\n", 
		    nifti_datatype_string(datatype));
    }
    
    return dest;
    
}