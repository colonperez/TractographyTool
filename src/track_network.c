/*
 *  track_network.c
 *  
 *
 *  Created by William Triplett on 5/1/10.
 * Modified by Luis M Colon-Perez 9/15/15
 */

#include "../../TrackTools_GUI/tt_current_version.h"
#include "track_network.h"
#include <string.h>
#include <math.h>

runtime_options opts;

int main (int argc, char **argv)
{
    nifti_image *nim = NULL;
    tv_header tv;
    FILE *fp_in = NULL, *fp_out = NULL;
    int ctr, i, tmpsize, seed_check;
    int *roi_data = NULL;
    int iobytes;
	
    
    /* the function that does the computation per fiber */
    int (*p_netwrk_fun)(nifti_image *, network *, track *) = NULL;
    
    /* network measures */
    network nw;
    track tr;
    
    opts.output_type = OUTPUT_AS_TAB_DELIMITED;
    opts.final_sort = SORT_BY_WEIGHTS;
	  opts.conn_method = 0;
    
    if (0 == get_options(argc, argv, &opts)) 
        exit(1);
    
    fprintf(stderr, "track_network (%d) starting...\n", kTT_CURRENT_VERSION);
    
    /*------------------------------------------------------------------*/
    /* Open the track file and make sure we can read it.                */
    
    memset(&tv, 0, sizeof(tv));
    
    fp_in = fopen(opts.tv_infile, "rb");
	
	
    if (fp_in == NULL) {
        fprintf(stderr, "Unable to open track file.\n");
        perror("Reason");
        return 1;
    }
    
    iobytes = fread(&tv, 1, sizeof(tv), fp_in);
    
    if (strcmp(tv.id_string, "TRACK") || tv.hdr_size != 1000) {
        fprintf(stderr, "Bad trackvis file: %s.\n", opts.tv_infile);
        return 1;
    } else if (iobytes != sizeof(tv)) {
        fprintf(stderr, "Unable to open output file: %s\n", opts.tv_infile);
        perror("Reason");
        return 1;
    }
    
    if (0 == tv.n_count) {
        tv.n_count = get_track_count(fp_in, &tv);
        if (tv.n_count == -1) {
            fprintf(stderr, "Unable to read %s, cannot determine number of tracks.\n", 
                    opts.tv_infile);
            fclose(fp_in);
            return 1;
        }
    }
    
    /*------------------------------------------------------------------*/
    /* open up the node file  */
    nim = nifti_image_read(opts.nii_filename, 1);
    
    /* something went wrong in nii library */
    if (nim == NULL) { return 1; }
	/*------------------------------------------------------------------*/
    /* Open the seeds file and make sure we can read it.               */
	
	
	char *seed_point_file = malloc(sizeof(char)*(strlen(opts.tv_infile)+60));
	sprintf(seed_point_file, "%s.seeds", opts.output_file);
	nw.fp_seeds_nt = fopen(seed_point_file, "w");
	
	if (nw.fp_seeds_nt == NULL) {
        fprintf(stderr, "Unable to open seeds file.\n");
        perror("Reason");
        return 1;
    }
	
	sprintf(seed_point_file, "%s.seeds", opts.tv_infile);
	nw.fp_seeds_tr = fopen(seed_point_file, "r");
	free(seed_point_file);
	
	if (nw.fp_seeds_tr == NULL) {
      seed_check = 0;
      opts.conn_method = 0;
  } else{ 
      seed_check = 1;
	}
	/*--------------------------------------------------------------------*/
	
	
    /* make sure dimensions of track file match node file
     and that the node file is not floating-point type */
    if (nim->dim[1] != tv.dim[0] || nim->dim[2] != tv.dim[1] || nim->dim[3] != tv.dim[2]) {
        fprintf(stderr, "Dimensions of ROI do not match track file:\n");
        fprintf(stderr, "ROI: [ %d, %d, %d ]; Track file: [ %d, %d, %d ]\n",
                nim->dim[1], nim->dim[2], nim->dim[3], 
                tv.dim[0], tv.dim[1], tv.dim[2]);
        return 1;
    } else if (opts.round_ok ==0 && (nim->datatype & (DT_FLOAT | DT_DOUBLE))) {
        fprintf(stderr, "Input ROI is of non-integral type.\n");
        fprintf(stderr, "Use the -r option to round ROI values.\n");
        return 1;
    }
    
    /* we are going to convert the roi data from whatever data type it is in
     to 4 byte integer. */
    roi_data = malloc(sizeof(int)*nim->dim[1]*nim->dim[2]*nim->dim[2]);
    if (roi_data == NULL) {
        fprintf(stderr, "Unable to allocate memory for roi data.\n");
        return 1;
    }
    
    if (! nii_recast_to_int32(nim, roi_data)) {
        fprintf(stderr, "Nifti datatype not supported\n");
        return 1;
    } else {
        free(nim->data);
        nim->data = (void *)roi_data;
        nim->datatype = DT_INT32;
        nim->nbyper = 4;
    }
    
    init_network(nim, &nw);
    
    /* allocating memory for the track data structure */
    tr.t_data = malloc(kTV_POINT_BUFSIZE*sizeof(float));
    if (tr.t_data == NULL) {
        fprintf(stderr, "Unable to allocate track data buffer.\n");
        return 1;
    }
    tr.t_data_conn = malloc(kTV_POINT_BUFSIZE*sizeof(float));
    if (tr.t_data_conn == NULL) {
        fprintf(stderr, "Unable to allocate memory for track portions.\n");
        return 1;
    }
    
    tr.tv_hdr = &tv;
    tr.write_track = &write_track_portion_nop;
    
    if (opts.term_or_pass == 1) {
      
        if (opts.save_matching == 1) {
            tr.write_track = &write_track_portion;
        }
        
        if (opts.conn_method == 0) {
            p_netwrk_fun = &node_intersect_passthrough;
        } else if (opts.conn_method == 1) {
            p_netwrk_fun = &node_intersect_luis_published;
        } else if (opts.conn_method == 2) {
            p_netwrk_fun = &node_intersect_cp;
        } else {
            fprintf(stderr, "Invalid conn_method. Cannot Continue. Please select H, C or Cp or type -h for help");
            exit(1);
        }
    } else {
        p_netwrk_fun = &node_intersect_terminal;
    }
    
    /*------------------------------------------------------------------*/
    /* read one track at a time from the file and look for terminal points
     that lie in a node region. */
    fprintf(stderr, "Track count in: %d\n", tv.n_count);
    ctr = 0;
    
    if (opts.term_or_pass == 1 && opts.save_matching == 1) {
        char *tmpbuf = malloc(sizeof(char)*strlen(opts.output_file)+4);
        if (tmpbuf == NULL) {
            fprintf(stderr, "Unable to allocate memory for tempfile name.\n");
            return 1;
        }
        
        sprintf(tmpbuf, "%s.trk", opts.output_file);
        
        fp_out = fopen(tmpbuf, "w+b");
        if (fp_out == NULL) {
            fprintf(stderr, "Unable to open track output file: %s.\n", tmpbuf);
            perror("Reason");
            return 1;
        }
        
        iobytes = fwrite(&tv, 1, sizeof(tv), fp_out);
        if (iobytes != sizeof(tv)) {
            fprintf(stderr, "Unable to write trackvis header.\n");
            perror("Reason");
            return 1;
        }
        tr.fp_out = fp_out;
        free(tmpbuf);
    }
    
    for (i=0; i<tv.n_count; i++) {
        int trksize;
        int nconns;
        
        iobytes = fread(&trksize, 1, sizeof(int), fp_in);
        
        if (iobytes != sizeof(int)) {
            fprintf(stderr, "Unable to read from track file.\n");
            perror("Reason");
            return 1;
        } else if (trksize*3 > kTV_POINT_BUFSIZE) {
            fprintf(stderr, "Track %d too large for buffer (track_size = %d), exiting.\n",
                    i, trksize);
            fclose(fp_in);
            return 1;
        }
        
        tr.n_points = trksize;
        tr.n_elements = trksize*3;
        
        tmpsize = (trksize * (3 + tv.n_scalars) + tv.n_properties) * sizeof(float);
        iobytes = fread(tr.t_data, 1, tmpsize, fp_in);
        if (iobytes != tmpsize) {
            fprintf(stderr, "Unable to read from track file.\n");
            perror("Reason");
            return 1;
        }
        
        track_interpolate_linear(nim, &tr);
        
        nconns = (*p_netwrk_fun)(nim, &nw, &tr);
        
        if (nconns > 0 && 1) {
            ctr += nconns;
        }
        
        if (i % 100000 == 0) {
            fprintf(stderr, "Track: %d of %d processed.\n", i, tv.n_count);
        }
    }
    
    output_network(&nw, opts.output_type);
    
    fclose(fp_in);
    fclose(nw.fp_seeds_tr);
	
    if (opts.save_matching == 1) {
        tv.n_count = ctr;
        fseek(fp_out, 0, SEEK_SET);
        fwrite(&tv, 1, sizeof(tv), fp_out);
        fclose(fp_out);
		fclose(nw.fp_seeds_nt);
		
		
    }
    
    nifti_image_free(nim);
    free_network(&nw);
    free(tr.t_data);
    if (tr.t_data_conn != NULL)
        free(tr.t_data_conn);
    //This is new, if user selects non-Hagmann method and seed file is not present then we stop the tracking with this error message
    fprintf(stderr, "Track count out: %d\n", ctr);
	if (seed_check ==1) {fprintf(stderr, "Seed file found, connectivity defined as specified by user\n");
		}else {
			fprintf(stderr, "Seed file WAS NOT found, connectivity method set to 'H' see usage for options\n");
		}

    
    return 0;
    
}

void track_interpolate_linear(nifti_image *nim, track *t)
{
    
    float *temp;
    float cp[3];
    float np[3];
    
    int n = t->n_points;
    int i=0, j=0;
    
    for (i=0; i<n; i++) {
        
        cp[0] = t->t_data[3*i+0]/nim->dx;
        cp[1] = t->t_data[3*i+1]/nim->dy;
        cp[2] = t->t_data[3*i+2]/nim->dz;
        
        t->t_data_conn[3*j+0] = t->t_data[3*i+0];
        t->t_data_conn[3*j+1] = t->t_data[3*i+1];
        t->t_data_conn[3*j+2] = t->t_data[3*i+2];
        ++j;
        
        if (i < (n-1)) {
            
            np[0] = t->t_data[3*(i+1)+0]/nim->dx;
            np[1] = t->t_data[3*(i+1)+1]/nim->dy;
            np[2] = t->t_data[3*(i+1)+2]/nim->dz;
            
            t->t_data_conn[3*j+0] = ((cp[0]+np[0])/2)*nim->dx;
            t->t_data_conn[3*j+1] = ((cp[1]+np[1])/2)*nim->dy;
            t->t_data_conn[3*j+2] = ((cp[2]+np[2])/2)*nim->dz;
            ++j;
            
        }
        
    }
    
    temp = t->t_data;
    t->t_data = t->t_data_conn;
    t->n_points = j;
    t->n_elements = j*3;
    t->t_data_conn = temp;
    memset(temp, 0, sizeof(float)*n*3);
    
}

int node_intersect_terminal(nifti_image *nim, network *nw, track *t)
{
    
    int conn_found = 0;
    int *roi_data = nim->data;
    int px, py, pz, index, node_1, node_2;
    int trksize = t->n_points;
    
    px = floor(t->t_data[3*0+0]/nim->dx);
    py = floor(t->t_data[3*0+1]/nim->dy);
    pz = floor(t->t_data[3*0+2]/nim->dz);
    
    index  = nii_voxel3_index(nim, 
                              FORCE_BETWEEN(px, 0, nim->dim[1]),
                              FORCE_BETWEEN(py, 0, nim->dim[2]),
                              FORCE_BETWEEN(pz, 0, nim->dim[3]));
    node_1 =  roi_data[index];
    
    if (node_1 != 0) {
        
        px = floor(t->t_data[3*(trksize-1)+0]/nim->dx);
        py = floor(t->t_data[3*(trksize-1)+1]/nim->dy);
        pz = floor(t->t_data[3*(trksize-1)+2]/nim->dz);
        
        index  = nii_voxel3_index(nim, 
                                  FORCE_BETWEEN(px, 0, nim->dim[1]),
                                  FORCE_BETWEEN(py, 0, nim->dim[2]),
                                  FORCE_BETWEEN(pz, 0, nim->dim[3]));
        node_2 = roi_data[index];
        
        if (node_2 != 0 && node_1 != node_2) {
            float len;
            nw->conn_count[node_1][node_2]++;
            nw->conn_count[node_2][node_1]++;
            len = get_track_length(t->tv_hdr, t->t_data, trksize);
            if (len != 0.0) {
                nw->weights[node_1][node_2] += 1.0/len;
                nw->weights[node_2][node_1] = nw->weights[node_1][node_2];
                nw->lengths[node_1][node_2] += len;
                nw->lengths[node_2][node_1] = nw->lengths[node_1][node_2];
            }
            ++conn_found;
        }
    }
    
    return conn_found;
    
}

int node_intersect_passthrough(nifti_image *nim, network *nw, track *t)
{
		//This is Hagmann's method
		int n = t->n_points;
		int i, index = 0;
		int conn_found = 0;
		int *roi_data = (int *)(nim->data);
		char edgeflag = 0;
		int in_roi = 0, last_roi = 0;
		float run_len = 0.0;
		float innode_len = 0.0;
		float step_len;
		float curr_pt[3];
		
		/* t->t_data is in units of voxel size (mm), so convert to voxels */
		curr_pt[0] = t->t_data[0]/nim->dx;
		curr_pt[1] = t->t_data[1]/nim->dy;
		curr_pt[2] = t->t_data[2]/nim->dz;
		
		/* look inside the ROI volume to see if the current point is on an 
		 active voxel (inside an ROI) */
		index = nii_voxel3_index(nim, 
								 floor(curr_pt[0]),
								 floor(curr_pt[1]),
								 floor(curr_pt[2]));
		
		
		/* roi[index] = the label number of the ROI, or 0 if outsize an ROI */
		if (roi_data[index] != 0) {
			in_roi = roi_data[index];
			edgeflag = 1;
			last_roi = roi_data[index];
		} else {
			in_roi = 0; /* roi_data[index]; */
			edgeflag = 0;
		}
		
		t->n_points_conn = 0;
		
		for (i=1; i<n; i++) {
			
			float next_pt[3];
			float step_vec[3];
			
			add_track_portion(t, i, t->n_points_conn);
			
			/* t->t_data is in units of voxel size (mm), so convert to voxels */
			next_pt[0] = t->t_data[3*i+0]/nim->dx;
			next_pt[1] = t->t_data[3*i+1]/nim->dy;
			next_pt[2] = t->t_data[3*i+2]/nim->dz;
			
			/* now, get the vector between the current and next point, but we
			 want it to be in units of voxel size (mm) not voxels, since we
			 will be computing track length. */
			step_vec[0] = (next_pt[0]-curr_pt[0])*nim->dx;
			step_vec[1] = (next_pt[1]-curr_pt[1])*nim->dy;
			step_vec[2] = (next_pt[2]-curr_pt[2])*nim->dz;
			
			/* compute the magnitude of the step vector */
			step_len = sqrt(step_vec[0]*step_vec[0]+
							step_vec[1]*step_vec[1]+
							step_vec[2]*step_vec[2]);
			run_len += step_len;
			
			curr_pt[0] = next_pt[0];
			curr_pt[1] = next_pt[1];
			curr_pt[2] = next_pt[2];
			
			index = nii_voxel3_index(nim, floor(curr_pt[0]),
									 floor(curr_pt[1]),floor(curr_pt[2]));
			
			if (roi_data[index] != in_roi) {	    
				
				if (in_roi == 0) {
					/* track was travelling in free space, now hit an ROI */
					
					if (edgeflag == 0) {
						/* track has encountered its first ROI */
						
						/* printf("Point %03d Entering ROI %d, [ %f, %f, %f ]\n", i, roi_data[index],
						 floor(p[0]), floor(p[1]), floor(p[2])); */
						
						in_roi   = roi_data[index];	
						edgeflag = 1;
						last_roi = roi_data[index];
						
					} else if (roi_data[index] != last_roi) {
						/* track has reached another ROI, thus creating an edge */
						
						int node_1 = roi_data[index]+1;
						int node_2 = last_roi+1;
						int ex = 0;
						
						/* printf("*** record_connection: A=%d, B=%d, dist=%f\n", 
						 last_roi, roi_data[index], run_len);  */
						/* printf("Point %03d Entering ROI %d, [ %f, %f, %f ]\n", i, roi_data[index],
						 floor(p[0]), floor(p[1]), floor(p[2])); */
						
						nw->conn_count[node_1-1][node_2-1]++;
						nw->conn_count[node_2-1][node_1-1]++;
						
						if (run_len != 0.0) {
							nw->weights[node_1-1][node_2-1] += 1.0/run_len;
							nw->weights[node_2-1][node_1-1] = nw->weights[node_1-1][node_2-1];
							nw->lengths[node_1-1][node_2-1] += run_len;
							nw->lengths[node_2-1][node_1-1] = nw->lengths[node_1-1][node_2-1];
						}
						
						add_track_portion(t, i, t->n_points_conn);
						t->n_points_conn++;	
						for (ex = 0; ex < opts.num_exclude_nodes; ex++) {
							if (nw->node_map[node_1-1] == opts.exclude_nodes[ex] ||
								nw->node_map[node_2-1] == opts.exclude_nodes[ex]) {
								/* This node is an exclude-nodes node */
								break;
							}
						}
						
						if (ex == opts.num_exclude_nodes) {
							/* Node was not found to be an exclude-nodes node,
							 we can write the track to the track file. */
							t->write_track(t);
							/* Luis: write_track is a pointer to either write_track_portion()
							 ( if -save-matching is enabled, or write_track_portion_nop()
							 if -save-matching is not enabled. )
							 */
							++conn_found;
						} else {
							t->n_points_conn = 0;
							t->n_points = 0;
						}
						
						innode_len = 0.0;
						edgeflag   = 1;
						in_roi     = roi_data[index];
						
					} else if (roi_data[index] == last_roi) {
						/* track as encountered the last roi again for some reason.
						 we will reset the length to zero and start the length
						 and saved track from zero (so that the part of the track
						 involved in the loop is discarded.) */
						run_len = 0;
						t->n_points_conn = 0;
					}
					
				} else if (roi_data[index] == 0) {
					
					/* track was within an ROI, now in free space */
					/* printf("Point %03d Leaving ROI %d, [ %f, %f, %f ]\n", i, in_roi,
					 floor(p[0]), floor(p[1]), floor(p[2])); */
					
					/* set edges to be allowed, set the last_roi to be the roi we're
					 leaving, reset the measurement length. */
					edgeflag   = 1;
					innode_len = 0;
					last_roi   = in_roi;
					in_roi     = 0;
					run_len    = 0.0;
					
					/* try to record the point that was last within the roi */
					add_track_portion(t, i-1, 0);
					t->n_points_conn = 1;
					add_track_portion(t, i, t->n_points_conn);
					t->n_points_conn++;	
					
				} else {
					
					/* this happens when nodes touch each other and the streamline
					 has moved from one node to the next without entering free
					 space */
					
					int node_1 = roi_data[index]+1;
					int node_2 = in_roi+1;
					int ex = 0;
					
					nw->conn_count[node_1-1][node_2-1]++;
					nw->conn_count[node_2-1][node_1-1]++;
					
					if (run_len != 0.0) {
						nw->weights[node_1-1][node_2-1] += 1.0/run_len;
						nw->weights[node_2-1][node_1-1] = nw->weights[node_1-1][node_2-1];
						nw->lengths[node_1-1][node_2-1] += run_len;
						nw->lengths[node_2-1][node_1-1] = nw->lengths[node_1-1][node_2-1];
					}
					
					/* record the last point and current point only */
					add_track_portion(t, i-1, 0);
					t->n_points_conn = 1;
					add_track_portion(t, i, t->n_points_conn);
					t->n_points_conn++;	
					
					for (ex = 0; ex < opts.num_exclude_nodes; ex++) {
						if (nw->node_map[node_1-1] == opts.exclude_nodes[ex] ||
							nw->node_map[node_2-1] == opts.exclude_nodes[ex]) {
							/* node is found to be an exclude-nodes node */
							break;
						}
					}
					
					if (ex == opts.num_exclude_nodes) {
						/* node is not an exclude-nodes node, write track */
						t->write_track(t);
						++conn_found;
					} else {
						t->n_points_conn = 0;
						t->n_points = 0;
					}
					
					innode_len = 0.0;
					edgeflag   = 1;
					in_roi     = roi_data[index];
					
				}
				
			} else if (in_roi != 0) {
				
				innode_len += step_len;
				t->n_points_conn++;
				
			} else if (in_roi == 0) {
				
				t->n_points_conn++;	
				
			}
			
		}
		
		return conn_found;
}

int node_intersect_luis_published(nifti_image *nim, network *nw, track *t)
{		
		//this is Colon-Perez method
		int n = t->n_points, m;
		int i, j,k=0, index = 0, flag_nw=0;
		int conn_found = 0;
		int *roi_data = (int *)(nim->data);
		char edgeflag = 0;
		int in_roi = 0, last_roi = 0;
		float trcO=0.0, trcT=0.0, trcTh=0.0 , discard=0 ;
		float run_len = 0.0;
		float innode_len = 0.0;
		float step_len;
		float curr_pt[3];
		float *ind=NULL;
		
		ind = malloc(n*3*sizeof(float));
        if (ind == NULL) {
			fprintf(stderr, "Unable to allocate memory for ind (size = %lu)\n",
					n*3*sizeof(float));
			exit(1);
        }
		fscanf(nw->fp_seeds_tr, "%f, %f , %f, %f \n", &trcO, &trcT, &trcTh, &discard);
		trcO = floor(trcO);
		trcT = floor(trcT);
		trcTh = floor(trcTh);
		/* t->t_data is in units of voxel size (mm), so convert to voxels */
		curr_pt[0] = t->t_data[0]/nim->dx;
		curr_pt[1] = t->t_data[1]/nim->dy;
		curr_pt[2] = t->t_data[2]/nim->dz;
		
		/* look inside the ROI volume to see if the current point is on an 
		 active voxel (inside an ROI) */
		index = nii_voxel3_index(nim, 
								 floor(curr_pt[0]),
								 floor(curr_pt[1]),
								 floor(curr_pt[2]));
		
		/* roi[index] = the label number of the ROI, or 0 if outsize an ROI */
		if (roi_data[index] != 0) {
			in_roi = roi_data[index];
			edgeflag = 1;
			last_roi = roi_data[index];
			flag_nw =1;
		} else {
			in_roi = 0; /* roi_data[index]; */
			edgeflag = 0;
			flag_nw =0;
		}
		
		
		t->n_points_conn = 0;
		
		for (i=1; i<n; i++) {
			
			float next_pt[3];
			float step_vec[3];
     
			add_track_portion(t, i, t->n_points_conn);
			
			/* t->t_data is in units of voxel size (mm), so convert to voxels */
			next_pt[0] = t->t_data[3*i+0]/nim->dx;
			next_pt[1] = t->t_data[3*i+1]/nim->dy;
			next_pt[2] = t->t_data[3*i+2]/nim->dz;
			
			/* now, get the vector between the current and next point, but we
			 want it to be in units of voxel size (mm) not voxels, since we
			 will be computing track length. */
			step_vec[0] = (next_pt[0]-curr_pt[0])*nim->dx;
			step_vec[1] = (next_pt[1]-curr_pt[1])*nim->dy;
			step_vec[2] = (next_pt[2]-curr_pt[2])*nim->dz;
			
			/* compute the magnitude of the step vector */
			step_len = sqrt(step_vec[0]*step_vec[0]+
							step_vec[1]*step_vec[1]+
							step_vec[2]*step_vec[2]);
			run_len += step_len;
			
			curr_pt[0] = next_pt[0];
			curr_pt[1] = next_pt[1];
			curr_pt[2] = next_pt[2];
			
			index = nii_voxel3_index(nim, floor(curr_pt[0]),
									 floor(curr_pt[1]),floor(curr_pt[2]));
			
			if (roi_data[index] != in_roi) {	    
				
				if (in_roi == 0) {
					/* track was travelling in free space, now hit an ROI */
					
					
					if (edgeflag == 0) {
						/* track has encountered its first ROI */
						/* printf("Point %03d Entering ROI %d, [ %f, %f, %f ]\n", i, roi_data[index],
						 floor(p[0]), floor(p[1]), floor(p[2])); */
						
						in_roi   = roi_data[index];	
						edgeflag = 1;
						last_roi = roi_data[index];
						//flag_nw =1;
						
						
					} else if (roi_data[index] != last_roi) {
						/* track has reached another ROI, thus creating an edge */
						int node_1 = roi_data[index]+1;
						int node_2 = last_roi+1;
						k = 0;
						int ex =0;
						/* printf("*** record_connection: A=%d, B=%d, dist=%f\n", 
						 last_roi, roi_data[index], run_len);  */
						/* printf("Point %03d Entering ROI %d, [ %f, %f, %f ]\n", i, roi_data[index],
						 floor(p[0]), floor(p[1]), floor(p[2])); */
						
						nw->conn_count[node_1-1][node_2-1]++;
						nw->conn_count[node_2-1][node_1-1]++;
						
						
						add_track_portion(t, i, t->n_points_conn);
						t->n_points_conn++;
						
						for (ex = 0; ex < opts.num_exclude_nodes; ex++) {
							if (nw->node_map[node_1-1] == opts.exclude_nodes[ex] ||
								nw->node_map[node_2-1] == opts.exclude_nodes[ex]) {
								/* This node is an exclude-nodes node */
								break;
							}
						}
						
						if (ex == opts.num_exclude_nodes) {
						/* since we have a complete track, and are about to check if we write it
						 compute mean and stdev angle */
						
						/* if I do this roi_data[trackPoint]==in_roi, then add streamline   */
						
						for (j=0; j<t->n_points_conn; j++) { 
							ind[3*j+0] = t->t_data_conn[3*j+0]/nim->dx;
							ind[3*j+1] = t->t_data_conn[3*j+1]/nim->dy;
							ind[3*j+2] = t->t_data_conn[3*j+2]/nim->dz;
						}
						flag_nw = 0;
						
						m= t->n_points_conn;
						
						for (k=1; k < ((t->n_points_conn)) ; k++) {
							if (trcO==floor(ind[3*k+0])) {
								if (trcT==floor(ind[3*k+1])) {
									if (trcTh==floor(ind[3*k+2])) {
										if (!(floor(ind[3*k+0])==floor(ind[0]) && floor(ind[3*k+1])==floor(ind[1]) && floor(ind[3*k+2])==floor(ind[2]))) {
											if (!(floor(ind[3*k+2])==floor(ind[3*m-1]) && floor(ind[3*k+1])==floor(ind[3*m-2]) && floor(ind[3*k+0])==floor(ind[3*m-3]))) {
												fprintf(nw->fp_seeds_nt, "%f, %f , %f \n", trcO, trcT, trcTh);	
												t->write_track(t);
												innode_len = 0.0;
												edgeflag   = 1;											
												++conn_found;
												flag_nw=1;
												break;
											}
										}
									}
								}
							}
						}
							
						} else {
							t->n_points_conn = 0;
							t->n_points = 0;
							//edgeflag = 0;
							//flag_nw = 0;
						}
						in_roi     = roi_data[index];
						if (flag_nw==1) {
							if (run_len != 0.0) {
								nw->weights[node_1-1][node_2-1] += 1.0/run_len;
								nw->weights[node_2-1][node_1-1] = nw->weights[node_1-1][node_2-1];
								nw->lengths[node_1-1][node_2-1] += run_len;
								nw->lengths[node_2-1][node_1-1] = nw->lengths[node_1-1][node_2-1];
							}
							
						}else {
							nw->conn_count[node_1-1][node_2-1]--;
							nw->conn_count[node_2-1][node_1-1]--;
							t->n_points_conn--;	
						}
						
					} else if (roi_data[index] == last_roi) {
						/* track as encountered the last roi again for some reason.
						 we will reset the length to zero and start the length
						 and saved track from zero (so that the part of the track
						 involved in the loop is discarded.) */
						run_len = 0;
						t->n_points_conn = 0;
					}
					
				} else if (roi_data[index] == 0) {
					
					/* track was within an ROI, now in free space */
					/* printf("Point %03d Leaving ROI %d, [ %f, %f, %f ]\n", i, in_roi,
					 floor(p[0]), floor(p[1]), floor(p[2])); */
					
					/* set edges to be allowed, set the last_roi to be the roi we're
					 leaving, reset the measurement length. */
					edgeflag   = 1;
					innode_len = 0;
					last_roi   = in_roi;
					in_roi     = 0;
					run_len    = 0.0;
					flag_nw    = 1;
					
					/* try to record the point that was last within the roi */
					add_track_portion(t, i-1, 0);
					t->n_points_conn = 1;
					add_track_portion(t, i, t->n_points_conn);
					t->n_points_conn++;
				} else {
					
					/* this happens when nodes touch each other and the streamline
					 has moved from one node to the next without entering free
					 space */
					
					int node_1 = roi_data[index]+1;
					int node_2 = in_roi+1;
					int ex =0;
					run_len    = 0.0;
					/* try to record the point that was last within the roi */
					add_track_portion(t, i-1, 0);
					t->n_points_conn = 1;
					
					nw->conn_count[node_1-1][node_2-1]++;
					nw->conn_count[node_2-1][node_1-1]++;
					
					
					add_track_portion(t, i, t->n_points_conn);
					t->n_points_conn++;	
					for (ex = 0; ex < opts.num_exclude_nodes; ex++) {
						if (nw->node_map[node_1-1] == opts.exclude_nodes[ex] ||
							nw->node_map[node_2-1] == opts.exclude_nodes[ex]) {
							/* This node is an exclude-nodes node */
							break;
						}
					}
					
					if (ex == opts.num_exclude_nodes) {
					for (j=0; j<t->n_points_conn; j++) { 
						ind[3*j+0] = t->t_data_conn[3*j+0]/nim->dx;
						ind[3*j+1] = t->t_data_conn[3*j+1]/nim->dy;
						ind[3*j+2] = t->t_data_conn[3*j+2]/nim->dz;
					}
					flag_nw =0;
					m= t->n_points_conn;
					for (k=1; k< ((t->n_points_conn)); k++) {
						if (trcO==floor(ind[3*k+0])) {
							if (trcT==floor(ind[3*k+1])) {
								if (trcTh==floor(ind[3*k+2])) {
									if (!(floor(ind[3*k+0])==floor(ind[0]) && floor(ind[3*k+1])==floor(ind[1]) && floor(ind[3*k+2])==floor(ind[2]))) {
										if (!(floor(ind[3*k+2])==floor(ind[3*m-1]) && floor(ind[3*k+1])==floor(ind[3*m-2]) && floor(ind[3*k+0])==floor(ind[3*m-3]))) {
											fprintf(nw->fp_seeds_nt, "%f, %f , %f \n", trcO, trcT, trcTh);	
											t->write_track(t);
											innode_len = 0.0;
											edgeflag   = 1;
											++conn_found;
											flag_nw=1;
											
											break;
										}
									}
								}
							}
						}
					}
				} else {
					t->n_points_conn = 0;
					t->n_points = 0;
					//edgeflag = 0;
					//flag_nw = 0;
				}
					in_roi     = roi_data[index];
					if (flag_nw==1) {
						if (run_len != 0.0) {
							nw->weights[node_1-1][node_2-1] += 1.0/run_len;
							nw->weights[node_2-1][node_1-1] = nw->weights[node_1-1][node_2-1];
							nw->lengths[node_1-1][node_2-1] += run_len;
							nw->lengths[node_2-1][node_1-1] = nw->lengths[node_1-1][node_2-1];
						}
						
						
					}else {
						nw->conn_count[node_1-1][node_2-1]--;
						nw->conn_count[node_2-1][node_1-1]--;
						t->n_points_conn--;	
					}
					
					
				}
				
			} else if (in_roi != 0) {
				
				innode_len += step_len;
				t->n_points_conn++;
				
			} else if (in_roi == 0) {
				
				t->n_points_conn++;	
			}
			
		}
		
		free(ind);
		
		return conn_found;

}

int node_intersect_cp(nifti_image *nim, network *nw, track *t)
{
		//this is new Cp method, where the streamlines within the ROI are included in the edge weight. 
	int n = t->n_points, m;
	int i, j,k=0, index = 0, flag_nw=0;
	int conn_found = 0;
	int *roi_data = (int *)(nim->data);
	char edgeflag = 0;
	int in_roi = 0, last_roi = 0;
	float trcO=0.0, trcT=0.0, trcTh=0.0 , discard=0 ;
	float run_len = 0.0;
	float innode_len = 0.0;
	float step_len;
	float curr_pt[3];
	float *ind=NULL;
	int index_seed =0;
	
	ind = malloc(n*3*sizeof(float));
	if (ind == NULL) {
		fprintf(stderr, "Unable to allocate memory for ind (size = %lu)\n",
				n*3*sizeof(float));
		exit(1);
	}
	fscanf(nw->fp_seeds_tr, "%f, %f , %f, %f \n", &trcO, &trcT, &trcTh, &discard);
	trcO = floor(trcO);
	trcT = floor(trcT);
	trcTh = floor(trcTh);
	/* t->t_data is in units of voxel size (mm), so convert to voxels */
	curr_pt[0] = t->t_data[0]/nim->dx;
	curr_pt[1] = t->t_data[1]/nim->dy;
	curr_pt[2] = t->t_data[2]/nim->dz;
	
	/* look inside the ROI volume to see if the current point is on an 
	 active voxel (inside an ROI) */
	index = nii_voxel3_index(nim, 
							 floor(curr_pt[0]),
							 floor(curr_pt[1]),
							 floor(curr_pt[2]));
	
	/* roi[index] = the label number of the ROI, or 0 if outsize an ROI */
	if (roi_data[index] != 0) {
		in_roi = roi_data[index];
		edgeflag = 1;
		last_roi = roi_data[index];
		flag_nw =1;
	} else {
		in_roi = 0; /* roi_data[index]; */
		edgeflag = 0;
		flag_nw =0;
	}
	
	
	t->n_points_conn = 0;
	
	for (i=1; i<n; i++) {
		
		float next_pt[3];
		float step_vec[3];
		
		add_track_portion(t, i, t->n_points_conn);
		
		/* t->t_data is in units of voxel size (mm), so convert to voxels */
		next_pt[0] = t->t_data[3*i+0]/nim->dx;
		next_pt[1] = t->t_data[3*i+1]/nim->dy;
		next_pt[2] = t->t_data[3*i+2]/nim->dz;
		
		/* now, get the vector between the current and next point, but we
		 want it to be in units of voxel size (mm) not voxels, since we
		 will be computing track length. */
		step_vec[0] = (next_pt[0]-curr_pt[0])*nim->dx;
		step_vec[1] = (next_pt[1]-curr_pt[1])*nim->dy;
		step_vec[2] = (next_pt[2]-curr_pt[2])*nim->dz;
		
		/* compute the magnitude of the step vector */
		step_len = sqrt(step_vec[0]*step_vec[0]+
						step_vec[1]*step_vec[1]+
						step_vec[2]*step_vec[2]);
		run_len += step_len;
		
		curr_pt[0] = next_pt[0];
		curr_pt[1] = next_pt[1];
		curr_pt[2] = next_pt[2];
		
		index = nii_voxel3_index(nim, floor(curr_pt[0]),
								 floor(curr_pt[1]),floor(curr_pt[2]));
		
		if (roi_data[index] != in_roi) {	    
			
			if (in_roi == 0) {
				/* track was travelling in free space, now hit an ROI */
				
				if (edgeflag == 0) {
					/* track has encountered its first ROI */
					/* printf("Point %03d Entering ROI %d, [ %f, %f, %f ]\n", i, roi_data[index],
					 floor(p[0]), floor(p[1]), floor(p[2])); */
					
					in_roi   = roi_data[index];	
					edgeflag = 1;
					last_roi = roi_data[index];
					//flag_nw =1;
					
					
				} else if (roi_data[index] != last_roi) {
					/* track has reached another ROI, thus creating an edge */
					int node_1 = roi_data[index]+1;
					int node_2 = last_roi+1;
					k = 0;
					int ex =0;
					/* printf("*** record_connection: A=%d, B=%d, dist=%f\n", 
					 last_roi, roi_data[index], run_len);  */
					/* printf("Point %03d Entering ROI %d, [ %f, %f, %f ]\n", i, roi_data[index],
					 floor(p[0]), floor(p[1]), floor(p[2])); */
					
					nw->conn_count[node_1-1][node_2-1]++;
					nw->conn_count[node_2-1][node_1-1]++;
					
					
					add_track_portion(t, i, t->n_points_conn);
					t->n_points_conn++;
					
					for (ex = 0; ex < opts.num_exclude_nodes; ex++) {
						if (nw->node_map[node_1-1] == opts.exclude_nodes[ex] ||
							nw->node_map[node_2-1] == opts.exclude_nodes[ex]) {
							/* This node is an exclude-nodes node */
							break;
						}
					}
					
					if (ex == opts.num_exclude_nodes) {
						
						for (j=0; j<t->n_points_conn; j++) { 
							ind[3*j+0] = t->t_data_conn[3*j+0]/nim->dx;
							ind[3*j+1] = t->t_data_conn[3*j+1]/nim->dy;
							ind[3*j+2] = t->t_data_conn[3*j+2]/nim->dz;
						}
						flag_nw = 0;
						
						m= t->n_points_conn;
						index_seed = nii_voxel3_index(nim, trcO, trcT, trcTh);
						if (roi_data[index_seed]==node_1-1||roi_data[index_seed]==node_2-1){
							fprintf(nw->fp_seeds_nt, "%f, %f , %f \n", trcO, trcT, trcTh);
							t->write_track(t);
							innode_len = 0.0;
							edgeflag   = 1;											
							++conn_found;
							flag_nw=1;
						}else {
						for (k=1; k < ((t->n_points_conn)) ; k++) {
							if (trcO==floor(ind[3*k+0])) {
								if (trcT==floor(ind[3*k+1])) {
									if (trcTh==floor(ind[3*k+2])) {
										if (!(floor(ind[3*k+0])==floor(ind[0]) && floor(ind[3*k+1])==floor(ind[1]) && floor(ind[3*k+2])==floor(ind[2]))) {
											if (!(floor(ind[3*k+2])==floor(ind[3*m-1]) && floor(ind[3*k+1])==floor(ind[3*m-2]) && floor(ind[3*k+0])==floor(ind[3*m-3]))) {
												fprintf(nw->fp_seeds_nt, "%f, %f , %f \n", trcO, trcT, trcTh);	
												t->write_track(t);
												innode_len = 0.0;
												edgeflag   = 1;											
												++conn_found;
												flag_nw=1;
												break;
								}	}	}	}	}	}	}	
					} else {
						t->n_points_conn = 0;
						t->n_points = 0;
					}
					in_roi     = roi_data[index];
					if (flag_nw==1) {
						if (run_len != 0.0) {
							nw->weights[node_1-1][node_2-1] += 1.0/run_len;
							nw->weights[node_2-1][node_1-1] = nw->weights[node_1-1][node_2-1];
							nw->lengths[node_1-1][node_2-1] += run_len;
							nw->lengths[node_2-1][node_1-1] = nw->lengths[node_1-1][node_2-1];
						}
						
					}else {
						nw->conn_count[node_1-1][node_2-1]--;
						nw->conn_count[node_2-1][node_1-1]--;
						t->n_points_conn--;	
					}
					
				} else if (roi_data[index] == last_roi) {
					/* track as encountered the last roi again for some reason.
					 we will reset the length to zero and start the length
					 and saved track from zero (so that the part of the track
					 involved in the loop is discarded.) */
					run_len = 0;
					t->n_points_conn = 0;
				}
				
			} else if (roi_data[index] == 0) {
				
				/* track was within an ROI, now in free space */
				/* printf("Point %03d Leaving ROI %d, [ %f, %f, %f ]\n", i, in_roi,
				 floor(p[0]), floor(p[1]), floor(p[2])); */
				
				/* set edges to be allowed, set the last_roi to be the roi we're
				 leaving, reset the measurement length. */
				edgeflag   = 1;
				innode_len = 0;
				last_roi   = in_roi;
				in_roi     = 0;
				run_len    = 0.0;
				flag_nw    = 1;
				
				/* try to record the point that was last within the roi */
				add_track_portion(t, i-1, 0);
				t->n_points_conn = 1;
				add_track_portion(t, i, t->n_points_conn);
				t->n_points_conn++;
			} else {
				
				/* this happens when nodes touch each other and the streamline
				 has moved from one node to the next without entering free
				 space */
				
				int node_1 = roi_data[index]+1;
				int node_2 = in_roi+1;
				int ex =0;
				run_len    = 0.0;
				/* try to record the point that was last within the roi */
				add_track_portion(t, i-1, 0);
				t->n_points_conn = 1;
				
				nw->conn_count[node_1-1][node_2-1]++;
				nw->conn_count[node_2-1][node_1-1]++;
				
				
				add_track_portion(t, i, t->n_points_conn);
				t->n_points_conn++;	
				for (ex = 0; ex < opts.num_exclude_nodes; ex++) {
					if (nw->node_map[node_1-1] == opts.exclude_nodes[ex] ||
						nw->node_map[node_2-1] == opts.exclude_nodes[ex]) {
						/* This node is an exclude-nodes node */
						break;
					}
				}
				
				if (ex == opts.num_exclude_nodes) {
					for (j=0; j<t->n_points_conn; j++) { 
						ind[3*j+0] = t->t_data_conn[3*j+0]/nim->dx;
						ind[3*j+1] = t->t_data_conn[3*j+1]/nim->dy;
						ind[3*j+2] = t->t_data_conn[3*j+2]/nim->dz;
					}
					flag_nw =0;
					m= t->n_points_conn;
					index_seed = nii_voxel3_index(nim, trcO, trcT, trcTh);
					if (roi_data[index_seed]==node_1-1||roi_data[index_seed]==node_2-1){
						fprintf(nw->fp_seeds_nt, "%f, %f , %f \n", trcO, trcT, trcTh);
						t->write_track(t);
						innode_len = 0.0;
						edgeflag   = 1;											
						++conn_found;
						flag_nw=1;
					}else {
					for (k=1; k< ((t->n_points_conn)); k++) {
						if (trcO==floor(ind[3*k+0])) {
							if (trcT==floor(ind[3*k+1])) {
								if (trcTh==floor(ind[3*k+2])) {
									if (!(floor(ind[3*k+0])==floor(ind[0]) && floor(ind[3*k+1])==floor(ind[1]) && floor(ind[3*k+2])==floor(ind[2]))) {
										if (!(floor(ind[3*k+2])==floor(ind[3*m-1]) && floor(ind[3*k+1])==floor(ind[3*m-2]) && floor(ind[3*k+0])==floor(ind[3*m-3]))) {
											fprintf(nw->fp_seeds_nt, "%f, %f , %f \n", trcO, trcT, trcTh);	
											t->write_track(t);
											innode_len = 0.0;
											edgeflag   = 1;
											++conn_found;
											flag_nw=1;
											break;
										}
									}
								}
							}
						}
					}
					}
				} else {
					t->n_points_conn = 0;
					t->n_points = 0;
					//edgeflag = 0;
					//flag_nw = 0;
				}
				in_roi     = roi_data[index];
				if (flag_nw==1) {
					if (run_len != 0.0) {
						nw->weights[node_1-1][node_2-1] += 1.0/run_len;
						nw->weights[node_2-1][node_1-1] = nw->weights[node_1-1][node_2-1];
						nw->lengths[node_1-1][node_2-1] += run_len;
						nw->lengths[node_2-1][node_1-1] = nw->lengths[node_1-1][node_2-1];
					}
					
					
				}else {
					nw->conn_count[node_1-1][node_2-1]--;
					nw->conn_count[node_2-1][node_1-1]--;
					t->n_points_conn--;	
				}
				
				
			}
			
		} else if (in_roi != 0) {
			
			innode_len += step_len;
			t->n_points_conn++;
			
		} else if (in_roi == 0) {
			
			t->n_points_conn++;	
		}
		
	}
	
	free(ind);
	
	return conn_found;	
}

int get_options (int argc, char **argv, runtime_options *opts)
{
    
    /*------------------------------------------------------------------*/
    /* parse the command line options */
    int optarg = 0;
    int i;
    
    if (argc < 2) {
        print_usage();
        exit(1);
    }
    
    for (i=1; i<argc; i++) {
        if (0 == strcmp(argv[i], "-r")) {
            opts->round_ok = 1;
            continue;
        } else if (0 == strcmp(argv[i], "-t")) {
            if ((i+1) != argc && 0 == strcmp(argv[i+1], "r_script")) {
                opts->output_type = OUTPUT_AS_R_SCRIPT;
            } else if ((i+1) != argc && 0 == strcmp(argv[i+1], "tab")) {
                opts->output_type = OUTPUT_AS_TAB_DELIMITED;
            } else {
                fprintf(stderr, "-t requires either 'r_script' or 'tab'.\n");
                return 0;
            }
            i++;
            continue;
        } else if (0 == strcmp(argv[i], "-s")) {
            if ((i+1) != argc && 0 == strcmp(argv[i+1], "weights")) {
                opts->final_sort = SORT_BY_WEIGHTS;
            } else if ((i+1) != argc && 0 == strcmp(argv[i+1], "lengths")) {
                opts->final_sort = SORT_BY_LENGTHS;
            } else if ((i+1) != argc && 0 == strcmp(argv[i+1], "num")) {
                opts->final_sort = SORT_BY_COUNT;
            } else {
                fprintf(stderr, "-s requires either 'weights', 'lengths', or 'num'.\n");
                return 0;
            }
            i++;
            continue;
			//This is new, user input to choose connectivity method
        } else if (0 == strcmp(argv[i], "-conn_method")) {
            if ((i+1) != argc && 0 == strcmp(argv[i+1], "H")) {
                opts->conn_method = 0;
            } else if ((i+1) != argc && 0 == strcmp(argv[i+1], "C")) {
                opts->conn_method = 1;
            } else if ((i+1) != argc && 0 == strcmp(argv[i+1], "Cp")) {
                opts->conn_method = 2;
            }else {
                fprintf(stderr, "-conn_method requires either 'H', 'C', or 'Cp': \n'H' for Hagmann, et al. Plos One(2007). DOI: 10.1371/journal.pone.0000597 \n'C' for Colon-Perez, et al. Plos One(2015). DOI: 10.1371/journal.pone.0131493 \n'Cp' similar to option 'C' but now including streamlines from within the ROI \n");
                return 0;
            }
            i++;
            continue;
		}else if (0 == strcmp(argv[i], "-m")) {
				if ((i+1) != argc && 0 == strcmp(argv[i+1], "pass")) {
					opts->term_or_pass = 1;
				} else if ((i+1) != argc && 0 == strcmp(argv[i+1], "term")) {
					opts->term_or_pass = 0;
				} else {
					fprintf(stderr, "-m requires either 'term' or 'pass'\n");
					return 0;
				}
				i++;
				continue;
        } else if (0 == strcmp(argv[i], "-save-matching")) {
            opts->save_matching = 1;
            continue;
        }	else if (0 == strcmp(argv[i], "-exclude-nodes")) {
            char *token = NULL, *arg_copy;
            const char delim[] = ",";
            int *exclude_nodes = malloc(sizeof(int)*25);
            int j = 0;
            
            memset(exclude_nodes, -1, sizeof(int)*25);
            
            arg_copy = strdup(argv[i+1]);
            token = strtok(arg_copy, delim);
            exclude_nodes[j++] = atoi(token); 
            while (token != NULL) {
                token = strtok(NULL, delim);
                if (token != NULL)
                    exclude_nodes[j++] = atoi(token);
            }
            opts->exclude_nodes = exclude_nodes;
            opts->num_exclude_nodes = j;
            i++;
            continue;    
        }
        if (optarg == 0) {
            opts->tv_infile = argv[i];
            optarg++;
        } else if (optarg == 1) {
            opts->nii_filename = argv[i];
            optarg++;
        } else if (optarg == 2) {
            opts->output_file = argv[i];
            optarg++;
        }
    }
    
    if (optarg != 3) {
        print_usage();
        exit(1);
    }    
    
    return 1;
}

/* convert x,y,z coordinate to single-dimension array index into nii array */
int nii_voxel3_index(nifti_image *nim,
                     const int x,
                     const int y,
                     const int z)
{
    return (x + y*nim->dim[1] + z*nim->dim[1]*nim->dim[2]);
}

void init_network(nifti_image *nim, network *nw)
{
    
    int i;
    
    /* count the number of regions in the roi file, so that each region becomes
     a node */
    fprintf(stderr, "Counting and preparing nodes... ");
    nw->node_count = count_nodes(nim, &(nw->node_map));
    if (nw->node_map == NULL) {
        fprintf(stderr, "unable to count nodes.\n");
        exit(1);
    }
    fprintf(stderr, "done.\n");
    
    /* allocate space for the connection matrix and the length, weights. */
    nw->conn_count = malloc(nw->node_count*sizeof(int*));
    nw->lengths    = malloc(nw->node_count*sizeof(float*));
    nw->weights    = malloc(nw->node_count*sizeof(float*));
    if (nw->conn_count == NULL || nw->lengths == NULL || nw->weights == NULL) {
        fprintf(stderr, "Unable to allocate memory for network accounting\n");
        fprintf(stderr, "in track_network.c:558.\n");
        exit(1);
    }
    for (i=0; i<nw->node_count; i++) {
        nw->conn_count[i] = malloc(nw->node_count*sizeof(int));
        nw->lengths[i] = malloc(nw->node_count*sizeof(float));
        nw->weights[i] = malloc(nw->node_count*sizeof(float));
	
        if (nw->weights[i] == NULL || nw->lengths[i] == NULL || nw->conn_count[i] == NULL) {
            fprintf(stderr, "Unable to allocate memory for network accounting\n");
            fprintf(stderr, "in track_network.c:567.\n");
            exit(1);
        }
        memset(nw->lengths[i], 0, nw->node_count*sizeof(float));
        memset(nw->conn_count[i], 0, nw->node_count*sizeof(int));
        memset(nw->weights[i], 0, nw->node_count*sizeof(float));
	
    }
}

void free_network(network *nw)
{
    
    int i;
    
    for (i=0; i<nw->node_count; i++) {
        free(nw->conn_count[i]);
        free(nw->lengths[i]);
        free(nw->weights[i]);
    }
    
    free(nw->conn_count);
    free(nw->lengths);
    free(nw->weights);
    free(nw->node_map);
    
}

/* conver nii file data to 4 byte int.*/
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

/* print out the usage information */
void print_usage (void)
{
	//This is new, slight modifications to usage output in terminal
    fprintf(stderr,
            "Usage: track_network [options] INPUT_TRACK_FILE ROI_NODE_FILE OUTPUT \n\n");
    fprintf(stderr, 
            "where [options] is given by:\n");
    fprintf(stderr,
            "   -t r_script    : output network data as an R script\n");
    fprintf(stderr,
            "   -t tab         : output network data as tab-delmited (default)\n");
    fprintf(stderr,
            "   -m term        : specify the method of edge calculation.\n");
    fprintf(stderr,
            "                    In order for a fiber to contribute to an edge, each\n");
    fprintf(stderr,
            "                    endpoint must lie within a node. (This is the default)\n");
    fprintf(stderr,
            "   -m pass        : A fiber may simply pass through a node to contribute to an.\n");
    fprintf(stderr,
            "                    edge. In this case, only the node-to-node length of the fiber\n");
    fprintf(stderr, 
            "                    is used to compute the edge weight.\n");    
    fprintf(stderr,
            "-exclude-nodes : set this option followed by a comma-separated (no spaces) list\n");
    fprintf(stderr,
            "                 of node numbers to exclude from trackvis file. The network results\n");
	fprintf(stderr,
			"                 will not be removed from the network csv file. The nodes will\n");
    fprintf(stderr,
            "                 be used in connectivity tracking, and will act as exclusion\n");
    fprintf(stderr,
            "                 nodes. Eg: -exclude-nodes 3,4,5\n");
	fprintf(stderr,
            "-save-matching    : Use this option in conjunction with '-m pass' to save only the\n");
    fprintf(stderr,
            "                   portions of the tracks that connect the nodes to a separate\n");
    fprintf(stderr,
            "                   file. The file will have the same name as OUTPUT with the\n");
    fprintf(stderr,
            "                   '.trk' extension appended.\n");
    fprintf(stderr,
            "   -s weights     : sort output by inverse lenght sum, weights corresponds to non-normalized edge weight (ascending) (default)\n");
    fprintf(stderr,
            "   -s num         : sort output by number of fibers (ascending)\n");
    fprintf(stderr,
            "   -s lengths     : sort output by average fiber length (ascending)\n\n");
	fprintf(stderr,
            "   -conn_method H : define connections as Hagmann,et.al.PlosOne(2007)DOI:10.1371/journal.pone.0000597(default)\n");
    fprintf(stderr,
            "   -conn_method C : define connections as Colon-Perez,et.al.Plos One(2015)DOI:10.1371/journal.pone.0131493\n");
    fprintf(stderr,
            "   -conn_method Cp: define connections similar to option 'C' but now including streamlines from within the ROI\n");
	fprintf(stderr,
            "   -conn_method    In order to use 'C' or 'Cp' which is a dimensionless and independent of resolution tracking measure, \n");
	fprintf(stderr,
            "                   the option -so in track_tracker has to be on. To do this just type -so as an option in track_traker \n");
    
    fprintf(stderr, 
            "track_network creates a graph network file from tracts in\n");
    fprintf(stderr,
            "INPUT_TRACK_FILE and nodes specified in ROI_NODE_FILE\n");
    fprintf(stderr,
            "(a nifti 3D volume). The network data is written to file OUTPUT\n");
    fprintf(stderr, 
            "which is created if it does not already exist.\n");
}

/* count the number of nodes in a nii node file.  */
int count_nodes (nifti_image *nim, int **node_map)
{
    
    int j, max, node_index, num_labels;
    int *bins, *reverse_bins;
    int dim = nim->nx*nim->ny*nim->nz;
    int *roi_data = nim->data;
    
    /* first find the max node number */
    max = 0;
    for (j=0; j<dim; j++) 
        if (roi_data[j] > max) 
            max = roi_data[j];
    
    /* we create <max> number of bins */
    bins         = malloc(sizeof(int)*(max+1));
    reverse_bins = malloc(sizeof(int)*(max+1));
    if (bins == NULL || reverse_bins == NULL) {
        fprintf(stderr, "Unable to allocate memory in count_nodes()\n");
        exit(1);
    }
    memset(bins, 0, sizeof(int)*(max+1));
    memset(reverse_bins, 0, sizeof(int)*(max+1));
    
    /* mark which bins are in use
     this ensures that we don't have problems
     with skipped label numbers in the roi file. starting at
     1 ensures that we ignore zero. */
    num_labels = 1;
    for (j=0; j<dim; j++) 
        if (roi_data[j] != 0 && bins[roi_data[j]]++ == 0)
            num_labels++;
    /* now num_labels contains the "true" number of bins, 
     excluding non-consecutively numbered labels */
    
    /* create a map from the integers to the roi node labels
     we will export this to the caller */
    *node_map = malloc(sizeof(int) * (num_labels));
    if (*node_map == NULL) {
        fprintf(stderr, "Unable to allocate memory for node map\n");
        exit(1);
    }
    
    (*node_map)[0] = 0;
    node_index = 1;
    for (j=1; j<=max; j++) {
        if (bins[j] != 0) {
            (*node_map)[node_index] = j;
            reverse_bins[j] = node_index;
            node_index++;
        }	    
    }
    
    /* now we replace the original segmentation labels with our node indices,
     which are numbered consecutively. The indices map into node_map and 
     identify the original segmentation number. For example:
     Suppose that in the segmentation file we have labels 1, 2, 10, and 12.
     Then at this point, node_map is a 5 element array, with 
     node_map[0] = 0; node_map[1] = 1; node_map[2] = 2;
     node_map[3] = 10; node_map[4] = 12.
     Now we will replace the original segmentation labels with the node
     map indices:
     */
    
    for (j=0; j<dim; j++) {
        if (roi_data[j] != 0) {
            roi_data[j] = reverse_bins[roi_data[j]];
        }
    }
    
    /* free up bins, but we leave node_map allocated since we return it as
     a side effect. */
    free(bins);
    free(reverse_bins);
    
    /* the number of nodes. */
    return node_index;
}

/* print the network info out as an R script to console. */
void output_network (network *nw, int type)
{
    
    int i,j,k;
    char *fmt_string;
    edge *edge_list;
    node *node_list;
    
    FILE *ofile;
    ofile = fopen(opts.output_file, "w");
	
    if (ofile == NULL) {
        fprintf(stderr, "Unable to open output file.\n");
        perror("Reason");
        exit(1);
    }        
    
    switch (type) {
        case OUTPUT_AS_R_SCRIPT:
            fprintf(ofile, "g <- network.initialize(%d, directed=FALSE);\n",
                    nw->node_count);
            fmt_string = "add.edge(g, %d, %d);  ## edge weight: %.4f, fib_count: %d, mean length: %.4f\n\0";
            break;
        case OUTPUT_AS_TAB_DELIMITED:
            fprintf(ofile, "Source Node\tDest Node\tInverse Sum of Streamline Lenght\tStreamline Count\tAverage Length\t\n");
            fmt_string = "%d\t%d\t%.4f\t%d\t%.4f\t\n\0";
            break;
        default:
            fclose(ofile);
            return;
            break;
    }
    
    edge_list = malloc(sizeof(edge)*(nw->node_count)*(nw->node_count+1)/2);
    if (edge_list == NULL) {
        fprintf(stderr, "Unable to allocate memory for edge list.\n");
        exit(1);
    }
    node_list = malloc(sizeof(node)*nw->node_count);
    if (node_list == NULL) {
        fprintf(stderr, "Unable to allocate memory for node list.\n");
        exit(1);
    } else {
        for (i=1; i<nw->node_count; i++) {
            node_list[i].degree = 0;
            node_list[i].strength = 0;
            node_list[i].surf_area = 0.0;
            node_list[i].node_num = nw->node_map[i];
        }
    }
    
    k=0;
    for (i=1; i<nw->node_count; i++) {
        
        for (j=1; j<i; j++) {
            
            if (nw->conn_count[i][j] != 0) {
                if (type == OUTPUT_AS_TAB_DELIMITED) {
                    edge_list[k].src_node = nw->node_map[i];
                    edge_list[k].dst_node = nw->node_map[j];
                } else {
                    edge_list[k].src_node = i;
                    edge_list[k].dst_node = j;
                }
                edge_list[k].weight = nw->weights[i][j];
                edge_list[k].n_fibers = nw->conn_count[i][j];
                edge_list[k].avg_fib_len = nw->lengths[i][j]/nw->conn_count[i][j];;
                
                node_list[i].degree++;
                node_list[i].strength += nw->weights[i][j];
                node_list[j].degree++;
                node_list[j].strength += nw->weights[i][j];
                
                k++;
            }
        }
    }
    
    switch (opts.final_sort) {
        case SORT_BY_WEIGHTS:
            qsort(edge_list, k, sizeof(edge), &edge_cmp_weights);
            break;
        case SORT_BY_COUNT:
            qsort(edge_list, k, sizeof(edge), &edge_cmp_count);
            break;
        case SORT_BY_LENGTHS:
            qsort(edge_list, k, sizeof(edge), &edge_cmp_lengths);
            break;
    }
    
    for (i=0; i<k; i++) {
        fprintf(ofile, fmt_string,
                edge_list[i].src_node, edge_list[i].dst_node,
                edge_list[i].weight, edge_list[i].n_fibers,
                edge_list[i].avg_fib_len);
    }
    
    if (type == OUTPUT_AS_TAB_DELIMITED) {
        fprintf(ofile, "\nNode\tDegree\tSurface Area\tStrength\n");
        for (i=1; i<nw->node_count; i++) {
            fprintf(ofile, "%04d\t%04d\t \t%.5f\n",
                    node_list[i].node_num,
                    node_list[i].degree,
                    node_list[i].strength);
        }
    }
    
    fclose(ofile);
    free(edge_list);
    free(node_list);
    
    return;
}

int edge_cmp_weights (const void *_e1, const void *_e2) {
    
    edge *e1 = (edge *) _e1;
    edge *e2 = (edge *) _e2;
    
    if (e1->weight == e2->weight) {
        return 0;
    } else {
        return e1->weight > e2->weight ? 1 : -1;
    }
    
}

int edge_cmp_lengths (const void *_e1, const void *_e2)
{
    
    edge *e1 = (edge *) _e1;
    edge *e2 = (edge *) _e2;
    
    if (e1->avg_fib_len == e2->avg_fib_len) {
        return 0;
    } else {
        return e1->avg_fib_len > e2->avg_fib_len ? 1 : -1;
    }
    
}

int edge_cmp_count (const void *_e1, const void *_e2) 
{
    
    edge *e1 = (edge *) _e1;
    edge *e2 = (edge *) _e2;
    
    if (e1->n_fibers == e2->n_fibers) {
        return 0;
    } else {
        return e1->n_fibers > e2->n_fibers ? 1 : -1;
    }
    
}

void add_track_portion (track *t, int src_idx, int dest_idx)
{
    if (t->t_data_conn == NULL) return;
    
	
    t->t_data_conn[3*dest_idx+0] = t->t_data[3*src_idx+0];
    t->t_data_conn[3*dest_idx+1] = t->t_data[3*src_idx+1];
    t->t_data_conn[3*dest_idx+2] = t->t_data[3*src_idx+2];
    
}    


/* Some macros to make things a little easier to follow */
#define VEC3_ASSIGN(LHS, RHS, INDEX) LHS[0] = RHS[3*(INDEX)+0]; \
				     LHS[1] = RHS[3*(INDEX)+1]; \
				     LHS[2] = RHS[3*(INDEX)+2];

#define VEC3_DOT_PRODUCT(VEC1, VEC2) (VEC1[0]*VEC2[0]+ \
				      VEC1[1]*VEC2[1]+ \
                                      VEC1[2]*VEC2[2])

/* *mean and *stdev are passed as pointers, so the function doesn't return
   anything it just modifies those two arguments */
void compute_angle_stats(track *t, float *mean, float *stdev)
{

    int i = 0;
    float last_pt[3];
    float curr_pt[3];
    float next_pt[3];
    
    float last_step_vec[3];
    float next_step_vec[3];

    /* t->t_data_conn holds the completed track that contributes to an edge. In
     * fact, t->t_data_conn has enough space for a track with 100,000 points. So,
     * we will find the index of the last point in the track and use the rest of
     * t->t_data_conn to hold the angles from which to compute the mean and stdev.
     */
    int angles_start = 3 * t->n_points_conn + 3;
    
    /* we will work from the center, so that the current point is the vertex
       of the angle made by the previous step and the next step */
    VEC3_ASSIGN(last_pt, t->t_data_conn, i)
    
    /* set these to 0, otherwise there will be trouble */
    *mean  = 0;
    *stdev = 0;
    
    for (i = 1; i < t->n_points_conn; i++) {
	
	/* make sure not to run past the last track point */
	if (i+1 == t->n_points_conn) break;
	
	VEC3_ASSIGN(curr_pt, t->t_data_conn, i)
	VEC3_ASSIGN(next_pt, t->t_data_conn, (i+1))
	
	last_step_vec[0] = last_pt[0]-curr_pt[0];
	last_step_vec[1] = last_pt[1]-curr_pt[1];
	last_step_vec[2] = last_pt[2]-curr_pt[2];
	
	next_step_vec[0] = curr_pt[0]-next_pt[0];
	next_step_vec[1] = curr_pt[1]-next_pt[1];
	next_step_vec[2] = curr_pt[2]-next_pt[2];
	
	float last_mag  = sqrt(VEC3_DOT_PRODUCT(last_step_vec, last_step_vec));
	float next_mag  = sqrt(VEC3_DOT_PRODUCT(next_step_vec, next_step_vec));
	float cos_theta = fabs(VEC3_DOT_PRODUCT(last_step_vec,next_step_vec)/(last_mag*next_mag));

	if (cos_theta > 1) { cos_theta = 1.0; }
	
	t->t_data_conn[angles_start+i] = acos(cos_theta)*180.0/M_PI;
	*mean += t->t_data_conn[angles_start+i];
	
    }
    
    /* here i = number of track points - 1 since there is one less angle than there
       are points */
    *mean /= (float)i;
    
    int j;
    for (j=1; j<=i ; j++) {
	*stdev += ((t->t_data_conn[angles_start+j] - *mean) * 
		   (t->t_data_conn[angles_start+j] - *mean));
    }
    *stdev = sqrt(*stdev/(float)i);

}


int write_track_portion (track *t)
{
    
    int nbytes;
    int bytes_out; 
    if (t == NULL) return 0;
    
    bytes_out = sizeof(float)*3*t->n_points_conn + sizeof(int);
    
    nbytes = fwrite(&(t->n_points_conn), 1, sizeof(int), t->fp_out);
    nbytes += fwrite(t->t_data_conn, 1, sizeof(float)*3*t->n_points_conn, t->fp_out);
    if (nbytes != bytes_out) {
        fprintf(stderr, "Unable to write track data to file.\n");
        perror("Reason");
        return 0;
    }
    
    t->n_points_conn = 0;
    
    return 1;
}

int write_track_portion_nop (track *t)
{
    return 1;
}
