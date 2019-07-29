/*
 *  track_track.h
 *  
 *
 *  Created by William Triplett on 5/9/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef TRACK_TRACK_H
#define TRACK_TRACK_H

#include "trackvis.h"
#include "nifti1_io.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef _WIN32
#define DIRECTORY_SEPARATOR '\\'
#else
#define DIRECTORY_SEPARATOR '/'
#endif

#define STREAM_DIRECTION_FW 2
#define STREAM_DIRECTION_BW 1
#define STREAM_DONE         0

#define STREAM_BUFSIZE      1000000
#define TRACK_BUFSIZE       STREAM_BUFSIZE*10
#define DEG2RAD(r)          r*180.0/M_PI

typedef enum nim_compare_result {
    EQUAL                 = 0,
    DIMENSION_MISMATCH    = 1,
    GRIDSPACE_MISMATCH    = 2,
    NUM_DIM_MISMATCH      = 4,
    ORIENT_QFORM_MISMATCH = 8,
    ORIENT_SFORM_MISMATCH = 16
} nim_compare_result;

typedef enum seed_plan_t {
    SEED_N_CUBED = 0,
    SEED_RANDOM  = 1
} seed_plan_t;

typedef struct streamline {
    int npts_fw;
    int npts_bw;
    int npts_ordered;
    int prop_direction;
    float seedpt[3];
    float currpt[3];
    float seeddir[3];
    float currdir[3];
    float *pts_fw;
    float *pts_bw;
    float *pts_ordered;
} streamline;

typedef struct track_params {
    nifti_image *seed_mask;
    nifti_image *term_mask;
    nifti_image *excl_mask;
    nifti_image *fa_image;
    nifti_image **direction_imgs;
    nifti_brick_list *direction_brks;
    int *loop_mask;
    char loop_check;
    char save_seeds;
    float **directions[5];
    int num_directions;
    float stream_len_thr;
    float fa_thr;
    float angle_thr;
    float cos_angle_thr;
    float voxel_dims[3];
    float step_size;
    float step_vec[3];
    int data_dims[3];
    int point_limit;
    int sd_dens;
    seed_plan_t seed_plan;
    char *track_buffer;
    int track_bufsize;
    int track_buf_used;
    char *seed_mask_file;
    char *term_mask_file;
    char *fa_mask_file;
    char *output_file;
    char *directions_dir;
    FILE *output_fp;
} track_params;

int track_perform_tracking(track_params *tp);

int nii_voxel3_index(nifti_image *nim,
		     const int x,
		     const int y,
		     const int z);

int nii_read_directions_from_files(char **file_names, track_params *tp);

nim_compare_result nii_compare_headers (nifti_image *nim1, nifti_image *nim2);

inline int get_direction_xyz(track_params *tp, int direction,
			     float *pt, float *result);
int get_direction_index(track_params *tp, int direction, int index,
			float *result);

void stream_reorganize_f2b(streamline *s);

int stream_check_fa_xyz(track_params *tp, float *pt);

int stream_check_term_xyz(track_params *tp, float *pt);

int stream_check_loop(track_params *tp, float *pt, int npts);

int stream_track(track_params *tp, streamline *s, int direction);
int stream_track_fw(track_params *tp, streamline *s);
int stream_track_bw(track_params *tp, streamline *s);

int nii_recast_to_int32 (nifti_image *nim);

float stream_find_best_direction(streamline *s, track_params *tp, float *best);

void track_flush_buffer(track_params *tp);

void track_add_stream_to_buffer(track_params *tp, streamline *s);

int track_sanity_check(track_params *tp);

void track_print_params(track_params *tp);

void track_print_usage(void);

#endif
