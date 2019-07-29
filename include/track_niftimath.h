/*
 *  track_niftimath.h
 *  track_tools
 *
 */

#ifndef TRACK_NIFTIMATH_H
#define TRACK_NIFTIMATH_H

#include "nifti1_io.h"

typedef enum _nm_operand_status {
    OPERAND_SCALAR          = 0,
    OPERAND_IMAGE_FOUND     = 1, 
    OPERAND_IMAGE_NOTFOUND  = 2,
    OPERAND_MISMATCH        = 4,
    OPERAND_UNKNOWN         = 8
} NM_OPERAND_STATUS;

typedef enum _nm_operation_types {
    NO_OP   = 0,
    OP_ADD  = 1,
    OP_SUB  = 2,
    OP_MUL  = 4,
    OP_DIV  = 8,
    OP_BIN  = 16,
    OP_THR  = 32,
    OP_UTHR = 64
} NM_OPERATION_TYPE;

typedef struct nifti_image_combo {
    nifti_image *nim;
    nifti_brick_list *nbl;
} nifti_image_combo;

typedef float comp_t;
const int comp_nifti_type = NIFTI_TYPE_FLOAT32;
#define NM_COMPUTATION_TYPE NIFTI_TYPE_FLOAT32;

void print_usage(void);

nifti_image_combo *nifti_read_combo(char *filename);

NM_OPERAND_STATUS get_operand_type (char *arg, NM_OPERATION_TYPE oper_type);
void nm_prepare_operand(nifti_image_combo *nic);
void nm_prepare_output(nifti_image_combo *nic, char *output_type);
void nm_free_operand(nifti_image_combo *nic);

int nii_recast_to_comp_t(nifti_image *nim, comp_t *data);

void nifti_combo_binary_op(nifti_image_combo *op1, 
                           nifti_image_combo *op2, 
                           NM_OPERATION_TYPE op_type);
void nifti_scalar_binary_op(nifti_image_combo *op1, 
                            float op2, 
                            NM_OPERATION_TYPE op_type);
void nifti_combo_unary_op(nifti_image_combo *op1,  
                          NM_OPERATION_TYPE op_type);


#endif /* TRACK_NIFTIMATH_H */

