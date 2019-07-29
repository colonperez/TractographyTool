/*
 *  track_density_mask.c
 *  
 *
 *  Created by William Triplett on 4/23/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "trackvis.h"
#include "nifti1_io.h"
#include "../../TrackTools_GUI/tt_current_version.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void print_usage (void);
int nii_recast_to_int32(nifti_image *nim, int *roi_data);
nifti_image *make_nii_from_reference_nii(char *ref_nim_filename);

int main (int argc, char **argv)
{
    nifti_image *dnim = NULL;
    tv_header *tvh = NULL;
    tv_file *tvf = NULL;
    tv_track *tvt = NULL;
    int ctr, i, j;
    int *roi_data = NULL;
    int px, py, pz, index, exitcode;
    char *tv_infile = NULL, *mask_filename = NULL;
    char *ref_nim_filename = NULL;
    
    int optarg = 0;
    int round_ok = 0;
    
    if (argc > 5 || argc < 3) {
        print_usage();
        return 1;
    }
    
    for (i=1; i<argc; i++) {
        if (0 == strcmp(argv[i], "-r")) {
            round_ok = 1;
            continue;
        }
        if (optarg == 0) {
            tv_infile = argv[i];
            optarg++;
        } else if (optarg == 1) {
            ref_nim_filename = argv[i];
            optarg++;
        } else if (optarg == 2) {
            mask_filename = argv[i];
            optarg++;
        }
    }
    
    if (optarg != 3) {
        print_usage();
        return 1;
    }
    
    fprintf(stderr, "track_density_mask (%d) starting...\n", kTT_CURRENT_VERSION);
    
    /* open input track file */
    tvf = tv_open(tv_infile);
    if (tvf == NULL) {
        return 1; 
    } else {
        tvh = tvf->tv_hdr;
    }
    
    fprintf(stderr, "Track count in: %d\n", tvh->n_count);
    
    /* allocate a new track data structure (trackvis.h) */
    tvt = tv_new_track();
    if (tvt == NULL) { 
        return 1;
    }
    
    dnim = make_nii_from_reference_nii(ref_nim_filename);
    if (dnim != NULL) {
        roi_data = (int *)dnim->data;
    } else {
        fprintf(stderr, "Unable to read reference NIFTI file.\n");
        return 1;
    }
    
    ctr = 0;
    exitcode = 0;
    
    for (i=0; i<tvh->n_count; i++) {
        
        int retcode = tv_read_next_track(tvf, tvt);
        if (retcode == 0) {
            goto CLEANUP;
        }
        
        int last_index = -1;
        
        for (j=0; j<tvt->n_points; j++) {
            
            px = floor(tvt->t_data[3*j+0]/tvh->voxel_size[0]);
            py = floor(tvt->t_data[3*j+1]/tvh->voxel_size[1]);
            pz = floor(tvt->t_data[3*j+2]/tvh->voxel_size[2]);
            /* px py pz are in voxel coordinates, tvh->voxel_order */
            
            if (px < 0) { px = 0; } else if (px >= dnim->dim[1]) { px = dnim->dim[1]-1; }
            if (py < 0) { pz = 0; } else if (py >= dnim->dim[2]) { py = dnim->dim[2]-1; }
            if (pz < 0) { pz = 0; } else if (pz >= dnim->dim[3]) { pz = dnim->dim[3]-1; }
            
            /* compute tvh->voxel_order to ROI voxel order here */
            
            index = (px + py*dnim->dim[1] + pz*dnim->dim[1]*dnim->dim[2]);
            
            if (index != last_index) {
                roi_data[index]++;
                last_index = index;
            }
        }
    }
    
CLEANUP:
    
    tv_close(tvf);
    if (mask_filename != NULL) {
        znzFile fp = znzopen(mask_filename, "wb", 0);
        free(dnim->fname);
        free(dnim->iname);
        dnim->fname = malloc(sizeof(char)*(1+strlen(mask_filename)));
        dnim->iname = malloc(sizeof(char)*(1+strlen(mask_filename)));
        if (dnim->fname == NULL || dnim->iname == NULL) {
            fprintf(stderr, "Unable to allocate memory for output file name.\n");
            return 1;
        }
        memcpy(dnim->fname, mask_filename, (1+strlen(mask_filename))*sizeof(char));
        memcpy(dnim->iname, mask_filename, (1+strlen(mask_filename))*sizeof(char));
        nifti_image_write_hdr_img2(dnim, 1, "wb", fp, NULL);
    }
    
    nifti_image_free(dnim);
    tv_free_track(tvt);
    
    return exitcode;
    
}

int nii_recast_to_int32 (nifti_image *nim, int *data)
{
    int i, nvoxels;
    
    nvoxels = nim->dim[1]*nim->dim[2]*nim->dim[3];
    
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
                fprintf(stderr, "Nifti datatype not supported: %s\n", 
                        nifti_datatype_string(nim->datatype));
                return 0;
        }
    }
    return 1;
}

nifti_image *make_nii_from_reference_nii(char *ref_nim_filename)
{
    nifti_image *nim, *out;
    
    /* open input ROI file and check dimensions are the same as track file */
    nim = nifti_image_read(ref_nim_filename, 0);
    if (nim == NULL) { 
        fprintf(stderr, "Unable to create nifti file.\n");
        return NULL;
    }
    
    out = nifti_copy_nim_info(nim);
    nifti_image_free(nim);
    nim = NULL;
    
    sprintf(out->descrip, "TrackTools Density Mask");
    
    out->datatype = NIFTI_TYPE_INT32;
    out->ndim     = 3;
    out->nbyper   = 4;
    out->nt       = 1;
    out->cal_max  = 0.0;
    out->cal_min  = 0.0;
    out->nvox     = out->nx*out->ny*out->nz;
    out->dim[4]   = 1;
    out->dim[0]   = 3;
    out->data     = malloc(4*out->nvox);
    if (out->data == NULL) {
        fprintf(stderr, "Unable to allocate memory for density mask.\n");
        return NULL;
    }
    
    memset(out->data, 0, 4*out->nvox);    
    return out;
}


void print_usage (void)
{
    fprintf(stderr,
            "Usage: track_density_mask INPUT_TRACK_FILE REFERENCE_NII OUTPUT_NII\n\n");
    fprintf(stderr, 
            "track_density_mask creates a nifti image where each voxel in\n");
    fprintf(stderr,
            "the image contains the number of tracks in INPUT_TRACK that\n");
    fprintf(stderr, 
            "pass through that voxel. The orientation information is taken\n");
    fprintf(stderr,
            "from REF_NII so that the dimensions and orientation are consistent\n");
    fprintf(stderr,
            "with the original data space. If OUTPUT_NII exists, it will be\n");
    fprintf(stderr,
            "overwritten silently.\n");
    
}

