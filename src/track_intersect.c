/*
 *  track_intersect.c
 *  
 *
 *  Created by William Triplett on 4/23/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "../../TrackTools_GUI/tt_current_version.h"
#include "trackvis.h"
#include "nifti1_io.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>

void print_usage (void);
int nii_recast_to_int32(nifti_image *nim, int *roi_data);

int main (int argc, char **argv)
{
    nifti_image *nim = NULL;
    tv_header *tvh = NULL;
    tv_file *tvf = NULL;
    tv_track *tvt = NULL;
    FILE *fp_out = NULL;
    FILE *fseed_pts_in = NULL;
    FILE *fseed_pts_out = NULL;
    int ctr, i, j, success;
    int *roi_data = NULL;
    int px, py, pz, index, iobytes, exitcode;
    char *tv_infile = NULL, *nii_filename = NULL, *tv_outfile = NULL;
    struct stat buf;
    
    int optarg = 0;
    int round_ok = 0;
    
    if (argc > 5 || argc < 3) {
        print_usage();
        return 1;
    }
    
    fprintf(stderr, "track_intersect (%d) starting...\n", kTT_CURRENT_VERSION);
    
    for (i=1; i<argc; i++) {
        if (0 == strcmp(argv[i], "-r")) {
            round_ok = 1;
            continue;
        }
        if (optarg == 0) {
            tv_infile = argv[i];
            optarg++;
        } else if (optarg == 1) {
            nii_filename = argv[i];
            optarg++;
        } else if (optarg == 2) {
            tv_outfile = argv[i];
            optarg++;
        }
    }
    
    if (optarg != 3) {
        print_usage();
        return 1;
    }
    
    if (0 == strcmp(tv_infile, tv_outfile)) {
        fprintf(stderr, "Input file and output track file cannot be the same.\n");
        return 1;
    }
    
    /* open input track file */
    tvf = tv_open(tv_infile);
    if (tvf == NULL) {
        return 1; 
    } else {
        tvh = tvf->tv_hdr;
    }
    
    /* Check to see if there is a seed point file, and if so open it for reading. */
    char *seed_point_file_in = malloc(sizeof(char)*(strlen(tv_infile)+10));
    sprintf(seed_point_file_in, "%s.seeds", tv_infile);
    if (stat(seed_point_file_in, &buf) == 0) {
	fseed_pts_in = fopen(seed_point_file_in, "r");
	fprintf(stderr, "Found seed point file: %s\n",seed_point_file_in);
    } else {
	fprintf(stderr, "No seed point file found (tried: %s)\n", seed_point_file_in);
    }
    free(seed_point_file_in);
    
    /* open input ROI file and check dimensions are the same as track file */
    nim = nifti_image_read(nii_filename, 1);
    if (nim == NULL) { return 1; }
    
    if (nim->dim[1] != tvh->dim[0] || 
        nim->dim[2] != tvh->dim[1] || 
        nim->dim[3] != tvh->dim[2]) {
        fprintf(stderr, "Dimensions of ROI do not match track file:\n");
        fprintf(stderr, "ROI: [ %d, %d, %d ]; Track file: [ %d, %d, %d ]\n",
                nim->dim[1], nim->dim[2], nim->dim[3], 
                tvh->dim[0], tvh->dim[1], tvh->dim[2]);
        return 1;
    } else if (round_ok ==0 && (nim->datatype & (DT_FLOAT | DT_DOUBLE))) {
        fprintf(stderr, "Input ROI is of non-integral type.\n");
        fprintf(stderr, "Use the -r option to round ROI values.\n");
        return 1;
    }
    
    /* we will recast whatever data type the ROI file is to 4byte int for ease */
    roi_data = malloc(sizeof(int)*nim->dim[1]*nim->dim[2]*nim->dim[2]);
    if (roi_data == NULL) {
        fprintf(stderr, "Unable to allocate memory for roi data.\n");
        return 1;
    } else {
        success = nii_recast_to_int32(nim, roi_data);
        if (! success) {
            fprintf(stderr, "Nifti datatype not supported\n");
            return 1;
        }
    }
    
    /* open the output track file */
    fp_out = fopen(tv_outfile, "w+b");
    if (NULL == fp_out) {
        fprintf(stderr, "Unable to open %s\n", tv_outfile);
        perror(NULL);
        return 1;
    }
    
    /* If there is a seed point input file, then open a seed point output file
       to write the seed points for tracks that intersect the ROI. */
    if (fseed_pts_in != NULL) {
	char *seed_point_file_out = malloc(sizeof(char)*(strlen(tv_outfile)+9));
	sprintf(seed_point_file_out, "%s.seeds", tv_outfile);
	fseed_pts_out = fopen(seed_point_file_out, "w");
	free(seed_point_file_out);
    }
    
    /* write the trackvis header, which will be similar to input tv header */
    iobytes = fwrite(tvh, 1, sizeof(*tvh), fp_out);
    if (iobytes != sizeof(*tvh)) {
        perror("Unable to write output file header: ");
        fclose(fp_out);
        return 1;
    }
    
    /* allocate a new track data structure (trackvis.h) */
    tvt = tv_new_track();
    if (tvt == NULL) { 
        return 1;
    }
    
    printf("Track count in: %d\n", tvh->n_count);
    ctr = 0;
    exitcode = 0;
    
    for (i=0; i<tvh->n_count; i++) {
        
	float seed_point[3];
	int tracksize;
	
        int retcode = tv_read_next_track(tvf, tvt);
	
        if (retcode == 0) {
	    fprintf(stderr, "Reading of track %d failed. Exiting.\n", i);
            goto CLEANUP;
        }
	
	/* Read a seed point for this track */
	if (fseed_pts_in != NULL)
	    fscanf(fseed_pts_in, "%f,%f,%f,%d\n", 
		   seed_point, seed_point+1, seed_point+2, &tracksize);
	
	
        for (j=0; j<tvt->n_points; j++) {
            
            px = floor(tvt->t_data[3*j+0]/tvh->voxel_size[0]);
            py = floor(tvt->t_data[3*j+1]/tvh->voxel_size[1]);
            pz = floor(tvt->t_data[3*j+2]/tvh->voxel_size[2]);
            /* px py pz are in voxel coordinates, tvh->voxel_order */
            
            if (px < 0) { px = 0; } else if (px >= nim->dim[1]) { px = nim->dim[1]-1; }
            if (py < 0) { pz = 0; } else if (py >= nim->dim[2]) { py = nim->dim[2]-1; }
            if (pz < 0) { pz = 0; } else if (pz >= nim->dim[3]) { pz = nim->dim[3]-1; }
            
            /* compute tvh->voxel_order to ROI voxel order here */
            
            index = (px + py*nim->dim[1] + pz*nim->dim[1]*nim->dim[2]);
            
            if ( roi_data[index] != 0 ) {
                iobytes = fwrite(&(tvt->n_points), 1, 4, fp_out);
                iobytes = fwrite(tvt->t_data, 1, tvt->n_elements*sizeof(float), fp_out);
		
		/* Write seed point for this track, if we have an open file. */
		if (fseed_pts_out != NULL) {
		    fprintf(fseed_pts_out, "%f,%f,%f,%d\n", 
			    seed_point[0], seed_point[1], seed_point[2], tracksize);
		}
		
                if (iobytes < tvt->n_elements) {
                    perror("Unable to write to output track file: ");
                    exitcode = 1;
                    goto CLEANUP;
                }
                
                ctr++;
                break;
            }		
            
        }
        
        if (i % 50000 == 0) {
            printf("Input tracks: %d, output tracks: %d, discarded: %d.\n",
                   i, ctr, i-ctr);
        }
        
    }
    
    tvh->n_count = ctr;
    fseek(fp_out, 0, SEEK_SET);
    iobytes = fwrite(tvh, 1, sizeof(*tvh), fp_out);
    if (iobytes != sizeof(*tvh)) {
        perror("Unable to write to output track file: ");
        exitcode = 1;
    }
    
    printf("Track count out: %d\n", ctr);
    
CLEANUP:
    fclose(fp_out);
    tv_close(tvf);
    nifti_image_free(nim);
    free(roi_data);
    tv_free_track(tvt);
    
    if (fseed_pts_out != NULL) {
	fclose(fseed_pts_out);
    }
    
    if (fseed_pts_in != NULL) {
	fclose(fseed_pts_in);
    }
    
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

void print_usage (void)
{
    fprintf(stderr,
            "Usage: track_intersect [-r] INPUT_TRACK_FILE ROI_FILE OUTPUT_TRACK_FILE\n\n");
    fprintf(stderr, 
            "tract_intersect searches INPUT_TRACT_FILE for tracts touching regions\n");
    fprintf(stderr,
            "specified in ROI_FILE (a nifti 3D volume) and saves matching tracts in\n");
    fprintf(stderr, 
            "OUTPUT_TRACT_FILE. If ROI_FILE contains floating point values, the -r\n");
    fprintf(stderr,
            "will optionally round them to integral values.\n");
}

