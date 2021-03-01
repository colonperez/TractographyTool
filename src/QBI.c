/*
 ============================================================================
Inputs:
e: m x 1 diffusion signal
Q: 3 x m column matrix of diffusion sampling wavevectors
	(e[i] is equal to the diffusion signal for the wavevector Q[i])
U: 3 x n matrix of reconstruction directions/directions of interest
V: 3 x p column matrix of basis function centers

LIST OF VARIABLES AND THEIR ASSOCIATED FUNCTIONS:
H: getDiffusionSignal
theta: getEquator (Not returned)
C: getEquator
S: getRotationMat
G: getDiffusionSignal
A: getReconstructionMatrix
Z: calculateNormalizationConstant
psi: computeODF
 ============================================================================
 */

#include"QBI.h"

#define DT_FLOAT32 16
#define NIFTI_FTYPE_NIFTI1_1  1

int main(int argc, char **argv) {
    QBI_RECON qbi;
    DIFF_DATA diff;

    qbi.diff = &diff;

    qbi_initialize_opts(&qbi, argc, argv);

    /* Set up the output data structure that will become the NIFTI direction
     files. */
    fprintf(stderr, "Initializing output data structures...\n");
    OUTPUT_DATA *output = initialize_output(diff.nii_image, qbi.num_output_files);

    int m = diff.n_b_high;
    int k = 30;

    /* Load the precomputed spherical domains */
    char *strbuf = malloc(sizeof(char)*strlen(qbi.datadir) + 15);

    fprintf(stderr, "Loading spherical domains...\n");
    sprintf(strbuf, "%s%c%s", qbi.datadir, DIRSEP, "tess_L1.dat");
    ICOS_TESS *restart_tess = load_tess_from_file(strbuf);
    
    sprintf(strbuf, "%s%c%s", qbi.datadir, DIRSEP, "tess_L3.dat");
    ICOS_TESS *reco_tess    = load_tess_from_file(strbuf);
    if (reco_tess == NULL || restart_tess == NULL) {
        fprintf(stderr, "Spherical tessellation files could not be loaded.\n");
        fprintf(stderr, "Make sure that the -datadir option points to the\n");
        fprintf(stderr, "directory that contains the 'tess_L*.dat' files.\n");
        exit(1);
    }
    free(strbuf); strbuf = NULL;
    
    qbi.reco_tess = reco_tess;  
    qbi.restart_tess = restart_tess;

    //Allocate Memory to Q, U, and V

    printf("Allocating memory...\n");

	double ** Q = malloc(sizeof(double*) * 3);
	for(int i = 0; i < 3; ++i){
		Q[i] = malloc(sizeof(double)*m);
	}

    int counter = 0;
	for(int i = 0; i < m; ++i){
        while(diff.bvals[counter]<200){
            counter++;
        }
        Q[0][i] = (double) diff.bvecs[counter];
        Q[1][i] = (double) diff.bvecs[diff.n_volumes+counter];
        Q[2][i] = (double) diff.bvecs[(2*diff.n_volumes)+counter];

        counter++;
	}

    int n = qbi.reco_tess->num_vertices;
    int n_reco_dirs = n;

    double** U = malloc(sizeof(double*)*3);
    double ** V = malloc(sizeof(double*) * 3);
    for(int i = 0; i < 3; ++i){
        U[i] = malloc(sizeof(double)*n);
        V[i] = malloc(sizeof(double)*n);
    }


    for(int i = 0; i < n; ++i){

        U[0][i] = qbi.reco_tess->vertices[i][0];
        V[0][i] = qbi.reco_tess->vertices[i][0];

        U[1][i] = qbi.reco_tess->vertices[i][1];
        V[1][i] = qbi.reco_tess->vertices[i][1];

        U[2][i] = qbi.reco_tess->vertices[i][2];
        V[2][i] = qbi.reco_tess->vertices[i][2];
    }

    double** A = malloc(sizeof(double*)*n);
    for(int i = 0; i < n; ++i){
        A[i] = malloc(sizeof(double)*m);
    }

    printf("Building Reconstruction Matrix...\n");

    //Inputs: 
        //Q : 3 × m column matrix of diffusion sampling wavevectors 
        //U : 3 × n column matrix of reconstruction points 
    //Returns the reconstruction matrix A, approximated using the "soft equator" A = ϕ(cos−1 |UT * Q|)
    getDiffusionSignal(U,Q,n,m,A);

    //Set up the maxima list
    MAXIMA *maxima_list  = malloc(sizeof(MAXIMA) * n_reco_dirs);
    if (maxima_list == NULL) {
        fprintf(stderr, "Unable to allocate memory for maxima list.\n");
        exit(1);
    }

	double* ODF = malloc(sizeof(double)*n);
    float* GFA_image = malloc(sizeof(float)*diff.nii_image->nz*diff.nii_image->ny*diff.nii_image->nx);

    fprintf(stderr, "Starting QBI reconstruction...\n");

    long unsigned int count = 0;
    for (int vz=0; vz<diff.nii_image->nz; vz++) {
        for (int vy=0; vy<diff.nii_image->ny; vy++) {
            for (int vx=0; vx<diff.nii_image->nx; vx++) {

                int n_maxima = 0;

                MAXIMA* max_list = malloc(sizeof(MAXIMA)*n);


            	int load_ok = load_voxel_double_highb(&diff, vx, vy, vz);

                if (-1 == load_ok) {
                    if (qbi.log_bad_voxels != 0) {
                        fprintf(stderr, "  WARNING: Voxel [%d,%d,%d] was not reconstructed (S0=0).\n",
                                vx, vy, vz);
                    }
                    continue;
                } else if (-2 == load_ok) {
                    if (qbi.log_bad_voxels != 0) {
                        fprintf(stderr, "  WARNING: Voxel [%d,%d,%d] was not reconstructed (nan/inf).\n",
                                vx, vy, vz);
                    }
                    continue;
                } else if (0 == load_ok) {
                    //Mask hit
                    count++;
                    continue;
                }

                double* e = diff.single_voxel_storage;

            	computeODF(A,e,n,m,ODF);

                n_maxima = find_local_maxima(qbi.reco_tess,ODF,qbi.prob_thresh,qbi.restart_tess,maxima_list);

                add_maxima_to_output(output, vx, vy, vz, U, maxima_list, n_maxima);

                //Get GFA value
                GFA_image[count] = (float) std(ODF,n)/rms(ODF,n);

                memset(ODF, 0, n*sizeof(double));

                free(max_list);
                count++;
                
            }
        }
        fprintf(stderr, "Slice: %d of %d Complete.\n", vz, diff.nii_image->nz);
        fflush(stderr);
    }

    if (qbi.GFAcompute == 1) {
        nifti_image *GFA_nim = nifti_simple_init_nim();
        memcpy(GFA_nim, diff.nii_image, sizeof(nifti_image));
        
        GFA_nim->datatype = DT_FLOAT32;
        GFA_nim->ndim     = 3;
        GFA_nim->nbyper   = 4;
        GFA_nim->nt       = 1;
        GFA_nim->nvox     = GFA_nim->nx * GFA_nim->ny * GFA_nim->nz;
        GFA_nim->dim[4]   = 1;
        GFA_nim->dim[0]   = 3;
        GFA_nim->data = GFA_image;
        GFA_nim->fname    = qbi.GFA_filename;
        GFA_nim->iname    = qbi.GFA_filename;
        GFA_nim->cal_max  = 0.0;
        GFA_nim->cal_min  = 0.0;
        qbi.GFA_filename = "GFA.nii.gz";

        fprintf(stderr, "Saving GFA image to %s.\n", qbi.GFA_filename);
        znzFile fp = znzopen(qbi.GFA_filename, "wb", 0);
        nifti_image_write_hdr_img2(GFA_nim, 1, "wb", fp, NULL);
    }

    fprintf(stderr, "QBI Reconstruction complete... saving output...\n");
    save_output(qbi.output_directory, output);

    fprintf(stderr, "Done.\n");

	for(int i = 0; i < 3; ++i){
		free(Q[i]);
		free(U[i]);
		free(V[i]);
	}
	free(U);
	free(Q);
	free(V);

    for(int i = 0; i < n; ++i){
        free(A[i]);
    }
    free(A);

    free(ODF);
    free(GFA_image);

	return 0;
}


void getA(double ** Q, double ** U, double ** V, int k, int m, int n, int p, double ** A) {

    double ** C = malloc(sizeof(double*)*3);
    for(int i = 0; i < 3; ++i){
        C[i] = malloc(sizeof(double)*k);
    }

    double** S = malloc(sizeof(double* )*3);
    for(int i = 0; i < 3;++i){
        S[i] = malloc(sizeof(double)*k*n);
    }

    double** G = malloc(sizeof(double*)*k*n);
    for(int i = 0; i < (k*n); ++i){
        G[i] = malloc(sizeof(double)*p);
    }

    double ** H = malloc(sizeof(double*) * m);
    for(int i = 0; i < m; ++i){
        H[i] = malloc(sizeof(double)*p);
    }

    getDiffusionSignal(Q,V,m,p,H); //Outputs H (ConvolutionAndPrediction)

    getEquator(k,C); //Outputs C (DirOfInterest)

    int counter = 0;
    double normalDirection[3] = {0,0,1};
    for(int i = 0; i < n; ++i){
        //For each direction of interest, multiply C by the rotation matrix and add it to S
        double Ui[3] = {U[0][i],U[1][i],U[2][i]};

        double** rotationMatrix = malloc(sizeof(double*)*3);
        for(int _ = 0; _ < 3; ++_){
            rotationMatrix[_] = malloc(sizeof(double)*3);
        }

        if(Ui[0]==0&&Ui[1]==0&&Ui[2]==-1){
            Ui[2] = 1;
        }

        getRotationMat(normalDirection,Ui,rotationMatrix);

        for(int j = 0; j < 3; ++j){
            for(int z = 0; z < k; z++){
                double sum = 0;
                for(int q = 0; q < 3; ++q){
                    sum += rotationMatrix[j][q]*C[q][z];
                }
                S[j][z+counter] = sum;
            }
        }

        counter += k;


        for(int _ = 0; _ < 3; ++_){
            free(rotationMatrix[_]);
        }
        free(rotationMatrix);

    }


    getDiffusionSignal(S,V,k*n,p,G); //Outputs G (ConvolutionAndPrediction)

    getReconstructionMatrix(G,H,k,n,p,m,A); //Outputs A (ConvolutionAndPrediction)

    for(int i = 0; i < 3; ++i){
        free(C[i]);
        free(S[i]);
    }
    free(C);
    free(S);

    for(int i = 0; i < (k*n); ++i){
        free(G[i]);
    }
    free(G);

    for(int i = 0; i < m; ++i){
        free(H[i]);
    }
    free(H);
}

void qbi_initialize_opts(QBI_RECON *qbi, int argc, char **argv)
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
    qbi->deco_p           = (double) DECO_P;
    qbi->deco_evals[0]    = (double) DECO_EVALS[0];
    qbi->deco_evals[1]    = (double) DECO_EVALS[1];
    qbi->deco_evals[2]    = (double) DECO_EVALS[2];
    qbi->prob_thresh      = (double) PROB_THRESH;
    qbi->num_output_files = 5;
    qbi->log_bad_voxels   = 0;
    qbi->datadir          = ".";
    qbi->data_filename    = NULL;
    qbi->bvec_filename    = NULL;
    qbi->mask_filename    = NULL;

    for (opt=1; opt<argc; opt++) {
        if (0 == strcmp(argv[opt], "-h")) {
            qbi_print_usage();
            exit(0);
        }
        if (0 == strcmp(argv[opt], "-odir")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -odir requires an argument.\n");
                exit(1);
            }
            qbi->output_directory = argv[opt+1];
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-log-bad-voxels")) {
            qbi->log_bad_voxels = 1;
            continue;
        } else if (0 == strcmp(argv[opt], "-pthresh")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -pthresh requires an argument.\n");
                exit(1);
            }
            qbi->prob_thresh = fabs(atof(argv[opt+1]));
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-mask")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -mask requires an argument.\n");
                exit(1);
            }
            qbi->mask_filename = argv[opt+1];
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-datadir")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -datadir requires an argument.\n");
                exit(1);
            }
            qbi->datadir = argv[opt+1];
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-ndir")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -ndir requires an argument.\n");
                exit(1);
            }
            qbi->num_output_files = atoi(argv[opt+1]);
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-data")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -data requires an argument.\n");
                exit(1);
            }
            qbi->data_filename = argv[opt+1];
            //printf(qbi->data_filename);
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-bvec")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -bvec requires an argument.\n");
                exit(1);
            }
            qbi->bvec_filename = argv[opt+1];
            opt++;
            continue;
        } 
        else if (0 == strcmp(argv[opt], "-S0")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -S0 requires an argument.\n");
                exit(1);
            }
            qbi->S0_filename = argv[opt+1];
            opt++;
            continue;
        }
        else if (0 == strcmp(argv[opt], "-GFAcompute")) {
            qbi->GFAcompute = 1;
            continue;
        }
        else if (0 == strcmp(argv[opt], "-bval")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -bval requires an argument.\n");
                exit(1);
            }       
            qbi->bval_filename = argv[opt+1];
            opt++;
            continue;
        }
        else {
            fprintf(stderr, "Ignoring junk on command line: %s\n", argv[opt]);
            fprintf(stderr, "Use -h to see a list of options.\n");
            exit(1);
        }
    }

    if (qbi->data_filename == NULL) {
        fprintf(stderr, "The -data <path> option is required to specify the location\n");
        fprintf(stderr, "of the diffusion-weighted data.\n");
        exit(1);
    } else if (qbi->bvec_filename == NULL) {
        fprintf(stderr, "The -bvec <path> option is required to specify the location\n");
        fprintf(stderr, "of the gradient directions for this data set.\n");
        exit(1);
    }
    if (qbi->S0_filename == NULL){
        printf("The -S0 <path> option is required.");
        exit(1);
    }

    fprintf(stderr, "S0 image will be loaded from %s.\n", qbi->S0_filename);
    qbi->diff->S0 = nifti_image_read(qbi->S0_filename, 1);
    if (qbi->diff->S0 == NULL) {
        exit(1);
    }

    /* these functions call exit() if anything bad happens. */
    fprintf(stderr, "Loading diffusion data...\n");
    read_diff_data_from_file(qbi->data_filename, qbi->diff);
    fprintf(stderr, "Loading gradient vectors...\n");
    read_bvecs_from_file(    qbi->bvec_filename, qbi->diff);
    fprintf(stderr, "Loading bvalues...\n");
    read_bvals_from_file(    qbi->bval_filename, qbi->diff);

    fprintf(stderr, "Loading binary mask...\n");
    qbi->diff->mask = nifti_image_read(qbi->mask_filename, 1);
    if (qbi->diff->mask == NULL) {
        fprintf(stderr, "Unable to load mask file: %s\n", qbi->mask_filename);
        fprintf(stderr, "I will attempt to reconstruct ALL voxels, even those\n");
        fprintf(stderr, "that may lie outside of the subject matter. Use results\n");
        fprintf(stderr, "with caution.\n\n");
        //exit(1);
    } else {
        nii_recast_to_int32(qbi->diff->mask);
    }

}

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

int read_acsii_file_to_float_array(char *filename, float *data, int data_size){

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
    //printf("1\n");
    free(strbuf);
    //printf("2\n");
    free(tmp);
    //printf("3\n");

    return elements_read;
}

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

///////////////////////////////////////////// Use gradient ascent to compute maximas on the sphere
int grad_find_max(double* ODF, int n, MAXIMA* maxima_list, int adj[199][4]){

    int numStartingPoints = 20;
    int startingIndices[numStartingPoints];
    for(int i = 0; i < numStartingPoints; ++i){ //Evenly sample starting points along the sphere
        startingIndices[i] = (i+1) * (n / numStartingPoints);
    }

    for(int i = 0; i < numStartingPoints; ++i){ //Perform gradient ascent at each starting point
        bool maxFound = false;
        int prevIndex = startingIndices[i];
        double prevVal = ODF[startingIndices[i]];

        //printf("NEW STARTING POINT \n\n");

        while(maxFound == false){ //Traverse the sphere until a maxima is found

            int dir1 = adj[prevIndex][0];
            int dir2 = adj[prevIndex][1];
            int dir3 = adj[prevIndex][2];
            int dir4 = adj[prevIndex][3]; //Indexes for each of the four directions

            double val1 = ODF[dir1];
            double val2 = ODF[dir2];
            double val3 = ODF[dir3];
            double val4 = ODF[dir4]; //Values for each of the four directions

            //printf("current: %d adjacent: %d(%f), %d(%f), %d(%f), %d(%f)\n",prevIndex,dir1,val1,dir2,val2,dir3,val3,dir4,val4);

            if(prevVal >= val1 && prevVal >= val2 && prevVal >= val3 && prevVal >= val4){ //If the previous value is greater than all of the traversal directions, a maxima has been found.
                maxFound = true;
                break;
            }
            else{ //If not, traverse in the direction that has the maximum increase
                if(val1 >= val2 && val1 >= val3 && val1 >= val4){
                    prevVal = val1;
                    prevIndex = dir1;
                }
                else if(val2 >= val1 && val2 >= val3 && val2 >= val4){
                    prevVal = val2;
                    prevIndex = dir2;
                }
                else if(val3 >= val1 && val3 >= val2 && val3 >= val4){
                    prevVal = val3;
                    prevIndex = dir3;
                }
                else if(val4 >= val1 && val4 >= val2 && val4 >= val3){
                    prevVal = val4;
                    prevIndex = dir4;
                }
                else{
                }
            }

        }

        //Once a maxima is found, add it to the list.

        maxima_list[i].index = prevIndex;
        maxima_list[i].value = prevVal;

    }

    //Filter out any redundant maxima

    int n_maxima = numStartingPoints; //Change based on if any maxima are filtered out

    for(int i = 0; i < n_maxima; ++i){ //Remove duplicate indexes
        for(int j = i+1; j < n_maxima; ++j){
            if(maxima_list[i].index == maxima_list[j].index){
                for(int k = j; k < n_maxima; ++k){
                    maxima_list[k] = maxima_list[k+1];
                }
                n_maxima--;
                j--;
            }
        }
    }

    qsort((void *)maxima_list, n_maxima, sizeof(MAXIMA), &maxima_compare);

    printf("First 5 Maxima values: \n");
    for(int i = 0; i < 5; ++i){
        printf("%f at %d\t",maxima_list[i].value,maxima_list[i].index);
    }
    printf("\n");

    return n_maxima;
}

int nii_voxel3_index(nifti_image *nim, const int x, const int y, const int z)
{
    return (x + y*nim->nx + z*nim->nx*nim->ny);
}

void qbi_print_usage(void)
{

    fprintf(stderr, "Usage: qbi_recon [OPTIONS]\n\n");
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
    fprintf(stderr, "  -ndir <integer>   Specifies the number of directions to save\n");
    fprintf(stderr, "                    per voxel. Maxima will be sorted by probability\n");
    fprintf(stderr, "                    and saved from highest probability to lowest.\n\n");
    fprintf(stderr, "  -pthresh <float>  Specifies the minimum probability (from 0 to 1)\n");
    fprintf(stderr, "                    that a maximal direction must have in order to be saved.\n\n");
    fprintf(stderr, "  -radius <float>   The displacement radius to consider.\n\n");
    fprintf(stderr, "  -odir <dirpath>   Path to a directory where the output files\n");
    fprintf(stderr, "                    should be saved. They will be named 'V_xx_all.nii\n");
    fprintf(stderr, "                    where xx ranges from 0 to ndir-1.\n\n");
    fprintf(stderr, "  -datadir <path>   path to the TrackTools data directory.\n\n");
    fprintf(stderr, "  -log-bad-voxels   If this option is used, qbi_recon will print\n");
    fprintf(stderr, "                    a list of all voxels that could not be reconstructed.\n\n");
    fprintf(stderr, "  -S0   <path>      path to an image containing the S0 data.\n");
    fprintf(stderr, "                    This is optional, and if not supplied the\n");
    fprintf(stderr, "                    S0 data will be computed internallly.\n\n");
    exit(0);

}

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
        sprintf(out[i]->descrip, "TrackTools QBI Reconstruction");

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
 * Adds the directions specified in maxima_list to the output data structure
 * at voxel location (x,y,z). Note, each maxima coordinate is added to a brick
 * in the brick list. That is, (x -> brick 0, y -> brick 1, z-> brick 2) and
 * the process is repeated for each maxima.
 *****************************************************************************/
void add_maxima_to_output(OUTPUT_DATA *output, int x, int y, int z, double **vertlist, MAXIMA *maxima_list, int n_maxima)
{
    int i;
    int index = nii_voxel3_index(output->nim[0], x, y, z);

    if (n_maxima > output->num_images) {
        /* fprintf(stderr, "  WARNING: Voxel: [%d,%d,%d], adding more maxima than number of files (%d > %d).\n",
         x, y, z, n_maxima, output->num_images); */
        n_maxima = output->num_images;
    }

    for (i=0; i<n_maxima; i++) {
        //float *mvert = vertlist[ maxima_list[i].index ];
        //((float *)output->nbl[i]->bricks[0])[index] = mvert[0];
        //((float *)output->nbl[i]->bricks[1])[index] = mvert[1];
        //((float *)output->nbl[i]->bricks[2])[index] = mvert[2];

        ((float *)output->nbl[i]->bricks[0])[index] = vertlist[0][maxima_list[i].index];
        ((float *)output->nbl[i]->bricks[1])[index] = vertlist[1][maxima_list[i].index];
        ((float *)output->nbl[i]->bricks[2])[index] = vertlist[2][maxima_list[i].index];
    }

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
	    fprintf(stderr, "Nifti datatype not supported.\n");
    }

    return dest;

}

int nii_recast_to_int32 (nifti_image *nim)
{
    int i, nvoxels;
    int *data;

    nvoxels = nim->dim[1]*nim->dim[2]*nim->dim[3];
    data = malloc(sizeof(int)*nvoxels);
    if (data == NULL) {
        fprintf(stderr, "Unable to allocate memory for seed mask.\n");//NOTE: Removed reference to nifti_datatype_string to run on my computer, can probably add back once I transfer the files to linux
        exit(1);
    }

    for (i=0; i<nvoxels; i++) {
        switch (nim->datatype) {
            case DT_UINT8:
                data[i] = (int) ((unsigned char *)(nim->data))[i];
                break;
            case DT_INT8:
                data[i] = (int) ((char *)(nim->data))[i];
                break;
            case DT_UINT16:
                data[i] = (int) ((unsigned short *)(nim->data))[i];
                break;
            case DT_INT16:
                data[i] = (int) ((short *)(nim->data))[i];
                break;
            case DT_UINT32:
                data[i] = (int) ((unsigned int *)(nim->data))[i];
                break;
            case DT_INT32:
                data[i] = (int) ((int *)(nim->data))[i];
                break;
            case DT_FLOAT32:
                data[i] = (int) round(((float *)(nim->data))[i]);
                break;
            case DT_FLOAT64:
                data[i] = (int) round(((double *)(nim->data))[i]);
                break;
            default:
                fprintf(stderr, "Nifti datatype not supported.\n"); //NOTE: Removed reference to nifti_datatype_string to run on my computer, can probably add back once I transfer the files to linux
                return 0;
        }
    }

    free(nim->data);
    nim->data = (void *)data;

    return 1;
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
