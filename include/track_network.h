/*
 *  track_network.h
 *  
 *
 *  Created by William Triplett on 5/7/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef TRACK_NETWORK_H
#define TRACK_NETWORK_H

#include "trackvis.h"
#include "nifti1_io.h"
#include <stdlib.h>
#include <math.h>

#define FORCE_BETWEEN(var, min, max) \
((var) < min) ? min : ( ((var) >= max) ? max - 1 : (var) )

#define OUTPUT_AS_R_SCRIPT 1
#define OUTPUT_AS_TAB_DELIMITED 2

#define SORT_BY_WEIGHTS 1
#define SORT_BY_LENGTHS 2
#define SORT_BY_COUNT   3

typedef struct track {
    float *t_data;      /* contains the track data */
    float *t_data_conn; /* contains the track data that has been "chopped" to connected nodes */
    int n_points_conn;  /* contains the number of (x,y,z) point in the chopped track */
    int n_points;       /* contains the number of (x,y,z) points */
    int n_elements;     /* contains the number of floats. generally 3*n_points */
    tv_header *tv_hdr;  /* contains the TrackVis header from whence the track came. */
    FILE *fp_out;
    int (*write_track)(struct track *);
} track;

typedef struct edgelist {
    void **edges;       /*       */
    int n_edges;        /*       */
} edgelist;

typedef struct edge {
    int src_node;       /* one node that defines this edge */
    int dst_node;       /* the other node that defines the edge */
    int n_fibers;       /* the number of fibers that contributed to this edge */
    float weight;       /* the Hagmann edge weight */
    float avg_fib_len;  /* the average length of the fibers contributing to edge */
} edge;

typedef struct node {
    int node_num;
    int degree;
    float strength;
    float surf_area;
    edgelist *el;
} node;

typedef struct runtime_options {
    int output_type;    /* OUTPUT_AS_R_SCRIPT or OUTPUT_AS_TAB_DELIMITED */
    int term_or_pass;   /* 0 = analyze using terminal points;
			   1 = analyze using passthrough method */
	int conn_method;   /* 0 = use Hagmman definition of connectivity(default);
						 1 = use Colon-Perez def of connectivity 
						2 = use Colon-Perez plus inside ROI def of connectivity*/
    int save_matching;   /* 0 = track portions that contribute to a network edge
			       will be NOT saved.
			   1 = track portions that contribute to a network edge
			       will be saved. */
    int round_ok;       /* OK to round floating point node .nii.gz files */
    int final_sort;     /* the field to sort the final output */
    char *tv_infile;    /* the name of the input trackvis file */
    char *exclude_nodes;
	int num_exclude_nodes;    /* number of the above */
    char *nii_filename; /* the name of the input nifti node file */
    char *output_file;  /* the name of the output text file */
} runtime_options;

typedef struct network {
    int *node_map;      /* a single dim array mapping indices to segmentation 
			   label numbers in the original nii node file */
    int node_count;     /* the number of distinct nodes */
    float **lengths;    /* 2 dim array containing the sum of the lengths of 
			   fibers connecting node [a][b] */
    float **weights;    /* 2 dim array conaining the Hagmann edge weights of
			   fibers connecting node [a][b] */
    int **conn_count;   /* 2 dim array containing the number of fibers connecting
			   node [a][b] */
	FILE *fp_seeds_tr;
	 FILE *fp_seeds_nt;
	
} network;

/*  prototypes */
void print_usage (void);
void track_interpolate_linear(nifti_image *nim, track *t);

int node_intersect_passthrough(nifti_image *nim, network *nw, track *t);
int node_intersect_terminal(nifti_image *nim, network *nw, track *t);
int node_intersect_luis_published(nifti_image *nim, network *nw, track *t);
int node_intersect_cp(nifti_image *nim, network *nw, track *t);

int nii_recast_to_int32(nifti_image *nim, int *roi_data);

void init_network(nifti_image *nim, network *nw);
void free_network(network *nw);

int count_nodes (nifti_image *nim, int **node_map);

void output_network (network *nw, int type);

int edge_cmp_weights(const void *e1, const void *e2);
int edge_cmp_lengths(const void *e1, const void *e2);
int edge_cmp_count(const void *e1, const void *e2);

int nii_voxel3_index(nifti_image *nim,
			    const int x,
			    const int y,
			    const int z);

int get_options(int argc, char **argv, runtime_options *opts);

void add_track_portion(track *t, int src_idx, int dest_idx);
int write_track_portion (track *t);
int write_track_portion_nop (track *t);



#endif
