/*
 *  track_niftimath.c
 *  track_tools
 *
 */

#include "../../TrackTools_GUI/tt_current_version.h"
#include "track_niftimath.h"
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/stat.h>

int main (int argc, char *argv[])
{
    
    int i, pass;
    NM_OPERAND_STATUS operand;
    NM_OPERATION_TYPE oper_type;
    nifti_image_combo *nic_in = NULL;
    
    for (i=0; i<argc; i++) {
        if (strcmp(argv[i], "-h") == 0) {
            print_usage();
            return 0;
        }
    }
    
    fprintf(stderr, "track_niftimath (%d) starting...\n", kTT_CURRENT_VERSION);
    
    for (pass=0; pass<2; pass++) {
        
        if (pass > 0) {
            nic_in = nifti_read_combo(argv[1]);
            if (nic_in == NULL) {
                fprintf(stderr, "Unable to read input nifti image.\n");
                return 1;
            } else {
                nm_prepare_operand(nic_in);
            }
        }
        
        for (i=2; i<argc; i++) {
            
            oper_type = NO_OP;
            
            /* find out which operation to perform */
            if        (strcmp(argv[i], "-add") == 0) {
                oper_type = OP_ADD;
            } else if (strcmp(argv[i], "-sub") == 0) {
                oper_type = OP_SUB;
            } else if (strcmp(argv[i], "-mul") == 0) {
                oper_type = OP_MUL;	    
            } else if (strcmp(argv[i], "-div") == 0) {
                oper_type = OP_DIV;
            } else if (strcmp(argv[i], "-thr") == 0) {
                oper_type = OP_THR;
            } else if (strcmp(argv[i], "-uthr") == 0) {
                oper_type = OP_UTHR;
            } else if (strcmp(argv[i], "-bin") == 0) {
                oper_type = OP_BIN;
                if (pass > 0) nifti_combo_unary_op(nic_in, oper_type);
                continue;
            } else {
                /* no operation, then it is time to save the output. */
                if (pass > 0) {
                    
                    if (i+2 < argc && strcmp(argv[i+1], "-odt") == 0)
                        nm_prepare_output(nic_in, argv[i+2]);
                    
                    sprintf(nic_in->nim->iname, "%s", argv[i]);
                    sprintf(nic_in->nim->fname, "%s", argv[i]);
                    znzFile fp = znzopen(nic_in->nim->iname, "wb", 0);
                    
                    if (fp != NULL) {
                        nifti_image_write_hdr_img2(nic_in->nim, 1, "wb", 
                                                   fp, nic_in->nbl);
                        return 0;
                    } else {
                        fprintf(stderr, "Unable to open output file.\n");
                        return 1;
                    }
                }
                break;
            }
            
            /* We made it here because oper_type != NO_OP, so we have
             * an operation to perform. Now, if there are no more arguments,
             * then something is wrong. */
            if (i+1 == argc) {
                fprintf(stderr, "Command line ended too early.\n");
                fprintf(stderr, "Expected operand after '%s'\n", argv[i]);
                return 1;
            }
            
            /* find out if the operand is valid for the operation, and
             * determine which type (scalar or image) the operand is. Then
             * perform the requested operation or error. */
            operand = get_operand_type(argv[i+1], oper_type);
            if (operand == OPERAND_SCALAR) {
                double op = strtod(argv[i+1], NULL);
                if (pass > 0)
                    nifti_scalar_binary_op(nic_in, (comp_t)op, oper_type);
                i++;
            } else if (operand == OPERAND_IMAGE_FOUND) {
                if (pass > 0) {
                    nifti_image_combo *nic_op = nifti_read_combo(argv[i+1]);
                    nm_prepare_operand(nic_op);
                    nifti_combo_binary_op(nic_in, nic_op, oper_type);
                    nm_free_operand(nic_op);
                }
                i++;
            } else if (operand == OPERAND_UNKNOWN) {
                fprintf(stderr, "Unknown operand for operator: %s\n", argv[i]);
                return 1;
            } else if (operand == OPERAND_MISMATCH) {
                fprintf(stderr, "Invalid operand for operator: %s\n", argv[i]);
                return 1;
            } else if (operand == OPERAND_IMAGE_NOTFOUND) {
                fprintf(stderr, "Cannot find image file for operand: %s\n",
                        argv[i]);
                return 1;
            }
        }
        
    }
    
    return 0;
	
}

NM_OPERAND_STATUS get_operand_type (char *arg, NM_OPERATION_TYPE oper_type)
{
    
    struct stat s;
    
    if (nifti_find_file_extension(arg) != NULL) {
        if (stat(arg, &s) == 0) {
            if (oper_type == OP_THR || oper_type == OP_UTHR) {
                return OPERAND_MISMATCH;
            } else {
                return OPERAND_IMAGE_FOUND;
            }
        } else {
            return OPERAND_IMAGE_NOTFOUND;
        }
    } else {
        int len = strlen(arg);
        int i;
        int seen_a_decimal = 0;
        int seen_a_sign    = 0;
        int seen_a_e       = 0;
        
        for (i=0; i<len; i++) {
            
            if (isdigit(arg[i])) continue;
            
            if (arg[i] == '+' || arg[i] == '-') {
                if (seen_a_sign) {
                    return OPERAND_UNKNOWN;
                } else {
                    seen_a_sign = 1;
                }
                continue;
            }
            
            if (arg[i] == '.') {
                if (seen_a_decimal) {
                    return OPERAND_UNKNOWN;
                } else {
                    seen_a_decimal = 1;
                }
                continue;
            }
            
            if (arg[i] == 'e' && i > 0 && i < (len-1)) {
                if (seen_a_e) {
                    return OPERAND_UNKNOWN;
                } else {
                    seen_a_e = 1;
                    if (arg[i+1] != '+' || arg[i+1] != '-') {
                        seen_a_sign = 1;
                        i++;
                    }
                    seen_a_decimal = 0;
                }
                continue;
            }
            
            return OPERAND_UNKNOWN;
            
        }
        
        return OPERAND_SCALAR;
    }
    
}

nifti_image_combo *nifti_read_combo(char *filename)
{
    
    nifti_brick_list *nbl_in = NULL;
    nifti_image *nim_in = NULL;
    nifti_image_combo *combo;
    
    combo = malloc(sizeof(nifti_image_combo));
    
    nbl_in = malloc(sizeof(nifti_brick_list));
    if (nbl_in == NULL) {
        fprintf(stderr, "Unable to allocate memory for nifti data.\n");
        return NULL;
    }
    
    nim_in = nifti_image_read_bricks(filename, 0, NULL, nbl_in);
    if (nim_in == NULL) { 
        fprintf(stderr, "Unable to read input data file.\n");
        return NULL;
        
    }
    
    combo->nim = nim_in;
    combo->nbl = nbl_in;
    
    return combo;
    
}

void nm_free_operand(nifti_image_combo *nic)
{
    
    int b;
    
    for (b=0; b<nic->nbl->nbricks; b++) {
        free(nic->nbl->bricks[b]);
    }
    
    nifti_image_free(nic->nim);
    free(nic->nbl);
}

void nm_prepare_operand(nifti_image_combo *nic)
{
    
    nifti_image *nim = nic->nim;
    nifti_brick_list *nbl = nic->nbl;
    int i;
    int nbricks = nbl->nbricks;
    int num_3d_vxl = nim->nx*nim->ny*nim->nz;
    comp_t *tmp;
    
    for (i=0; i<nbricks; i++) {
        
        nim->data = nbl->bricks[i];
        
        tmp = malloc(sizeof(comp_t) * num_3d_vxl);
        
        nii_recast_to_comp_t(nim, tmp);
        
        free(nbl->bricks[i]); nim->data = NULL;
        nbl->bricks[i] = tmp;
        
    }
    
    nim->datatype = comp_nifti_type;
    nim->nbyper   = sizeof(comp_t);
    nbl->bsize    = num_3d_vxl*sizeof(comp_t);
    
}

void nm_prepare_output(nifti_image_combo *nic, char *output_type)
{
    
    nifti_image *nim = nic->nim;
    nifti_brick_list *nbl = nic->nbl;
    int b,v;
    int nbricks = nbl->nbricks;
    int num_3d_vxl = nim->nx*nim->ny*nim->nz;
    
    void *tmp_data;
    
    int nifti_odt, nifti_odt_size;
    
    if (strcmp(output_type, "char") == 0) {
        nifti_odt = nifti_datatype_from_string("NIFTI_TYPE_INT8");
    } else if (strcmp(output_type, "short") == 0) {
        nifti_odt = nifti_datatype_from_string("NIFTI_TYPE_INT16");
    } else if (strcmp(output_type, "int") == 0) {
        nifti_odt = nifti_datatype_from_string("NIFTI_TYPE_INT32");
    } else if (strcmp(output_type, "float") == 0) {
        nifti_odt = nifti_datatype_from_string("NIFTI_TYPE_FLOAT32");
    } else if (strcmp(output_type, "double") == 0) {
        nifti_odt = nifti_datatype_from_string("NIFTI_TYPE_FLOAT64");
    } else {
        fprintf(stderr, "Output data type (%s) not supported.\n", output_type);
        fprintf(stderr, "Output will be saved with data type used for internal\n");
        fprintf(stderr, "calculations. Supported output types are:\n");
        fprintf(stderr, "    char, short, int, float, double\n");
        return;
    }
    
    nifti_datatype_sizes(nifti_odt , &nifti_odt_size, NULL);	
    
    for (b=0; b<nbricks; b++) {
        
        nim->data = nbl->bricks[b];
        
        tmp_data = malloc(nifti_odt_size*num_3d_vxl);
        
        for (v=0; v<num_3d_vxl; v++) {
            
            switch (nifti_odt) {
                    
                case DT_INT8:
                    (  (char*)tmp_data)[v] = (char)  ((comp_t*)nim->data)[v];
                    break;
                case DT_INT16:
                    ( (short*)tmp_data)[v] = (short) ((comp_t*)nim->data)[v];
                    break;
                case DT_INT32:
                    (   (int*)tmp_data)[v] = (int)   ((comp_t*)nim->data)[v];
                    break;
                case DT_FLOAT32:
                    ( (float*)tmp_data)[v] = (float) ((comp_t*)nim->data)[v];
                    break;
                case DT_FLOAT64:
                    ((double*)tmp_data)[v] = (double)((comp_t*)nim->data)[v];
                    break;
                default:
                    fprintf(stderr, "Output data type not supported.\n");
                    break;
                    
            }
        }
        
        nim->datatype = nifti_odt;
        nbl->bsize    = num_3d_vxl*nifti_odt_size;
        free(nbl->bricks[b]); nim->data = NULL;
        nbl->bricks[b] = tmp_data;
        
    }
    
}


void nifti_combo_binary_op(nifti_image_combo *op1, 
                           nifti_image_combo *op2, 
                           NM_OPERATION_TYPE oper_type)
{
    
    int b1, bconst, *b2;
    
    nifti_image *nim1 = op1->nim;
    nifti_brick_list *nbl1 = op1->nbl;
    nifti_brick_list *nbl2 = op2->nbl;
    
    int nvox = nim1->nx*nim1->ny*nim1->nz;
    
    if (nbl1->nbricks != nbl2->nbricks) {
        if (nbl2->nbricks != 1) {
            fprintf(stderr, "Number of volumes mismatch.\n");
            exit(1);
        } else {
            bconst = 0;
            b2 = &bconst;
        }
    } else {
        b2 = &b1;
    }
    
    for (b1 = 0; b1 < nbl1->nbricks; b1++) {
        int vox;
        for (vox=0; vox<nvox; vox++) {
            
            switch (oper_type) {
                case OP_ADD:
                    ((comp_t**)nbl1->bricks)[b1][vox] += 
                    ((comp_t**)nbl2->bricks)[*b2][vox];
                    break;
                case OP_SUB:
                    ((comp_t**)nbl1->bricks)[b1][vox] -= 
                    ((comp_t**)nbl2->bricks)[*b2][vox];
                    break;
                case OP_MUL:
                    ((comp_t**)nbl1->bricks)[b1][vox] *= 
                    ((comp_t**)nbl2->bricks)[*b2][vox];
                    break;
                case OP_DIV:
                    ((comp_t**)nbl1->bricks)[b1][vox] /= 
                    ((comp_t**)nbl2->bricks)[*b2][vox];
                    break;
                default:
                    break;
            }
        }
    }
}

void nifti_scalar_binary_op(nifti_image_combo *op1, 
                            comp_t op2, 
                            NM_OPERATION_TYPE oper_type)
{
    
    int b1;
    
    nifti_image *nim1 = op1->nim;
    nifti_brick_list *nbl1 = op1->nbl;
    
    int nvox = nim1->nx*nim1->ny*nim1->nz;
    
    for (b1 = 0; b1 < nbl1->nbricks; b1++) {
        int vox;
        for (vox=0; vox<nvox; vox++) {
            comp_t val;
            switch (oper_type) {
                case OP_ADD:
                    ((comp_t**)nbl1->bricks)[b1][vox] += op2;
                    break;
                case OP_SUB:
                    ((comp_t**)nbl1->bricks)[b1][vox] -= op2;
                    break;
                case OP_MUL:
                    ((comp_t**)nbl1->bricks)[b1][vox] += op2;
                    break;
                case OP_DIV:
                    ((comp_t**)nbl1->bricks)[b1][vox] /= op2;
                    break;
                case OP_THR:
                    val = ((comp_t**)nbl1->bricks)[b1][vox];
                    if (val < op2)
                        ((comp_t**)nbl1->bricks)[b1][vox] = 0;
                    break;
                case OP_UTHR:
                    val = ((comp_t**)nbl1->bricks)[b1][vox];
                    if (val > op2)
                        ((comp_t**)nbl1->bricks)[b1][vox] = 0;
                    break;
                default:
                    break;
            }
        }
    }
}

void nifti_combo_unary_op(nifti_image_combo *op1,  
                          NM_OPERATION_TYPE oper_type)
{
    
    int b1;
    
    nifti_image *nim1 = op1->nim;
    nifti_brick_list *nbl1 = op1->nbl;
    
    int nvox = nim1->nx*nim1->ny*nim1->nz;
    
    for (b1 = 0; b1 < nbl1->nbricks; b1++) {
        int vox;
        for (vox=0; vox<nvox; vox++) {
            comp_t val;
            switch (oper_type) {
                case OP_BIN:
                    val = ((comp_t**)nbl1->bricks)[b1][vox];
                    if (val != 0) 
                        ((comp_t**)nbl1->bricks)[b1][vox] = 1;
                    break;
                    
                default:
                    break;
            }
        }
    }
}

/* conver nii file data to internal type.*/
int nii_recast_to_comp_t (nifti_image *nim, comp_t *data)
{
    int i, nvoxels;
    
    nvoxels = nim->nx*nim->ny*nim->nz;
    
    for (i=0; i<nvoxels; i++) {
        switch (nim->datatype) {
            case DT_UINT8:
                data[i] = (comp_t) ((unsigned char *)(nim->data))[i];
                break;
            case DT_INT8:
                data[i] = (comp_t) ((char *)(nim->data))[i];
                break;
            case DT_UINT16:
                data[i] = (comp_t) ((unsigned short *)(nim->data))[i];
                break;
            case DT_INT16:
                data[i] = (comp_t) ((short *)(nim->data))[i];
                break;
            case DT_UINT32:
                data[i] = (comp_t) ((unsigned int *)(nim->data))[i];
                break;
            case DT_INT32:
                data[i] = (comp_t) ((int *)(nim->data))[i];
                break;
            case DT_FLOAT32:
                data[i] = (comp_t) ((float *)(nim->data))[i];
                break;
            case DT_FLOAT64:
                data[i] = (comp_t) ((double *)(nim->data))[i];
                break;
            default: 
                fprintf(stderr, "Nifti datatype not supported: %s\n", 
                        nifti_datatype_string(nim->datatype));
                return 0;
        }
    }
    return 1;
}

void print_usage(void)
{
    fprintf(stderr, "Usage: track_niftimath INPUT [[-OP [OPERAND]]...] OUTPUT [-odt OUTPUT_TYPE]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  where:\n");
    fprintf(stderr, "    INPUT is a nifti file, \n");
    fprintf(stderr, "    OP is an operation (see \"Operation Types\" below),\n");
    fprintf(stderr, "    OPERAND is either a nifti file, a scalar, or nothing depending on\n");
    fprintf(stderr, "            the type of OP,\n");
    fprintf(stderr, "    OUTPUT is the name of a nifti file that will hold the final output,\n");
    fprintf(stderr, "    OUTPUT_TYPE is the optional final data type of the output.\n");
    fprintf(stderr, "                (See \"Output Types\").\n\n");
    fprintf(stderr, "  Operation Types:\n");
    fprintf(stderr, "    -add  [file|scalar]  Voxelwise addition of file or scalar to INPUT\n");
    fprintf(stderr, "    -sub  [file|scalar]  Voxelwise addition of file or scalar to INPUT\n");
    fprintf(stderr, "    -mul  [file|scalar]  Voxelwise addition of file or scalar to INPUT\n");
    fprintf(stderr, "    -div  [file|scalar]  Voxelwise addition of file or scalar to INPUT\n");
    fprintf(stderr, "    -bin                 Binarize: replace nonzero voxels with 1\n");
    fprintf(stderr, "    -thr  [scalar]       Replace voxels less than scalar with 0.\n");
    fprintf(stderr, "    -uthr [scalar]       Replace voxels greater than scalar with 0.\n\n");
    fprintf(stderr, "    Note: operations may be chained together (see \"Examples\").\n");
    fprintf(stderr, "    Note: If INPUT has more than one volume, then file operands must\n");
    fprintf(stderr, "          have either one volume or the same number of volumes as INPUT.\n");
    fprintf(stderr, "          In the former case, the operation will be performed on each\n");
    fprintf(stderr, "          input volume. In the latter case, the operation will be performed\n");
    fprintf(stderr, "          volume-by-volume.\n\n");
    fprintf(stderr, "  Output Types:\n");
    fprintf(stderr, "    'char'   - 8-bit integer           'short' - 16-bit integer\n");
    fprintf(stderr, "    'int'    - 32-bit integer          'float' - 32-bit floating point\n");
    fprintf(stderr, "    'double' - 64-bit floating point.\n");
    fprintf(stderr, "    Note: all internal computations are done in 'float' type.\n\n");
    fprintf(stderr, "  Examples:\n");
    fprintf(stderr, "    track_niftimath input.nii.gz -add 3 -mul data1.nii.gz -bin output.nii.gz\n");
    fprintf(stderr, "       Add 3 to input.nii.gz, then multiply the result voxelwise by data1.nii.gz\n");
    fprintf(stderr, "       then binarize the result and save to output.nii.gz. Data will be saved\n");
    fprintf(stderr, "       in 32-bit floating point type. Note that there is no operand for -bin.\n\n");
    fprintf(stderr, "    track_niftimath input.nii.gz -thr 100 -uthr 200 output.nii.gz -odt char\n");
    fprintf(stderr, "       Take input.nii.gz and replace any voxels less than 100 with 0, then \n");
    fprintf(stderr, "       take the result of that and replace any voxels greater than 200 with 0\n");
    fprintf(stderr, "       and save the result to output.nii.gz in 8-bit integer format.\n\n");
    
}


