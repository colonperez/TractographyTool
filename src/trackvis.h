/*
 * trackvis.h
 *
 *  Created on: Mar 31, 2010
 *      Author: btt
 */

#ifndef TRACKVIS_H_
#define TRACKVIS_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include "nifti1_io.h"

#define kTV_POINT_BUFSIZE 300000

/*  See: http://www.trackvis.org/docs/?subsect=fileformat */
typedef struct tv_header {
    char id_string[6];
    short dim[3];
    float voxel_size[3];
    float origin[3];
    short n_scalars;
    char scalar_name[10][20];
    short n_properties;
    char property_name[10][20];
    float vox_to_ras[4][4];
    char reserved[444];
    char voxel_order[4];
    char pad2[4];
    float image_orientation_patient[6];
    char pad1[2];
    unsigned char invert_x;
    unsigned char invert_y;
    unsigned char invert_z;
    unsigned char swap_xy;
    unsigned char swap_yz;
    unsigned char swap_xz;
    int n_count;
    int version; 
    int hdr_size;
} tv_header;

typedef struct tv_file {
    tv_header *tv_hdr;
    FILE *tv_fp;
    int header_written;
    int read_only;
} tv_file;

typedef struct tv_track {
    float *t_data;      /* contains the track data */
    int n_points;       /* contains the number of (x,y,z) points */
    int n_elements;     /* contains the number of floats. generally 3*n_points */
} tv_track;

tv_file *tv_open(const char *filename);
void tv_close(tv_file *tvf);

tv_file *tv_create(const char *filename, tv_header *tvh);

tv_header *tv_clone_header(tv_header *tvh);

int tv_read_next_track(tv_file *tvf, tv_track *tvtr);

tv_track *tv_new_track(void);
void tv_free_track(tv_track *tvt);

int get_track_count(FILE *fp_in, tv_header *tv);
float get_track_length(tv_header *tv, float *track, int track_size);
void init_tv_header_from_nii(tv_header *tv, nifti_image *nii);

#endif /* TRACKVIS_H_ */
