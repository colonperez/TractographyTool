/*
 *  track_track.c
 *  
 *
 *  Created by William Triplett on 5/8/10.
 *
 */

#include "track_track.h"
#include "knuthrand.h"
#include <time.h>

int track_perform_tracking(track_params *tp)
{
    
    streamline stream;
    tv_header tv;
    float x, y, z;
    float *frac_seeds;
    int iobytes;
    FILE *fseed_pt = NULL;
    
    // zero out these data structures
    memset(&stream, 0, sizeof(stream));
    memset(&tv, 0, sizeof(tv));
    
    // set up the trackvis header using the seed mask.
    init_tv_header_from_nii(&tv, tp->seed_mask);
    
    // open up the trackvis file, write only
    tp->output_fp = fopen(tp->output_file, "w+b");
    if (tp->output_fp == NULL) {
        fprintf(stderr, "Unable to open output file: %s\n", tp->output_file);
        perror("Reason");
        exit(1);
    }
    
    if (tp->seed_plan == SEED_N_CUBED) {
        // create the fractional voxel seed points. These will be added to the 
        // voxel location to position the seed point within the voxel
        frac_seeds = malloc(sizeof(float)*tp->sd_dens);
        if (frac_seeds == NULL) {
            fprintf(stderr, "Unable to allocate memory for fractional seeds.\n");
            exit(1);
        }
        
        for (x=0; x<tp->sd_dens; x++) 
            frac_seeds[(int)x] = ((float)x + 0.5)/(float)tp->sd_dens;
        
    } else {
        long ctime = (long)time(NULL);
        fprintf(stderr, "Seeding randf with %ld\n", ctime);
        ranf_start(ctime);
    }
    
    // write out the trackvis header
    iobytes = fwrite(&tv, 1, sizeof(tv), tp->output_fp);
    if (iobytes != sizeof(tv)) {
        fprintf(stderr, "Unable to write to output file.\n");
        perror("Reason");
        exit(1);
    }
    
    fprintf(stderr, "Tracking...(no further output until finished.)\n");
    
    stream.pts_fw = malloc(sizeof(float) * STREAM_BUFSIZE);
    stream.pts_bw = malloc(sizeof(float) * STREAM_BUFSIZE);
    stream.pts_ordered = malloc(sizeof(float) * 2*STREAM_BUFSIZE);
    if (stream.pts_fw == NULL || stream.pts_bw == NULL || stream.pts_ordered == NULL) {
        fprintf(stderr, "Unable to allocate memory for streamline buffer.\n");
        exit(1);
    }
    
    /* Set up the seed point output file and open it for writing. It is just a
       plain text file. */
    if (tp->save_seeds) {
	char *seed_point_file = malloc(sizeof(char)*(strlen(tp->output_file)+10));
	sprintf(seed_point_file, "%s.seeds", tp->output_file);
	fseed_pt = fopen(seed_point_file, "w");
	free(seed_point_file);
    }
	
    // look for voxels that are "active" in the seed mask.
    for (z=0; z<tp->data_dims[2]; z++) {
        for (y=0; y<tp->data_dims[1]; y++) {
            for (x=0; x<tp->data_dims[0]; x++) {
                
                // this is the linear index into the seed mask data array
                int index = nii_voxel3_index(tp->seed_mask, x,y,z);
                
                if (((int *)tp->seed_mask->data)[index]) {
                    
                    // now track for each fractional seed point
                    int fx, fy, fz;
                    for(fx=0; fx<tp->sd_dens; fx++) {
                        for (fy=0; fy<tp->sd_dens; fy++) {
                            for (fz=0; fz<tp->sd_dens; fz++) {
                                int i;
                                
                                if (tp->seed_plan == SEED_N_CUBED) {
                                    stream.seedpt[0] = x + frac_seeds[fx];
                                    stream.seedpt[1] = y + frac_seeds[fy];
                                    stream.seedpt[2] = z + frac_seeds[fz];
                                } else {
                                    stream.seedpt[0] = x + (float)ranf_arr_next();
                                    stream.seedpt[1] = y + (float)ranf_arr_next();
                                    stream.seedpt[2] = z + (float)ranf_arr_next();
                                    fy = tp->sd_dens;
                                    fz = tp->sd_dens;
                                }
                                
                                
                                // track the multi-fiber system
                                for (i=0; i<tp->num_directions; i++) {
				    
                                    if (get_direction_index(tp, i, index,
                                                            stream.seeddir) == 0) {
                                        // direction is not valid (ie
                                        // zero vector)
                                        continue;
                                    }

                                    stream.npts_fw = 0;
                                    stream.npts_bw = 0;
                                    stream.prop_direction = STREAM_DIRECTION_FW;
                                    
                                    stream_track(tp, &stream, STREAM_DIRECTION_FW);
                                    stream_track(tp, &stream, STREAM_DIRECTION_BW);
                                    
                                    stream_reorganize_f2b(&stream);
                                    
				    /* Clean the loop check mask for the next
				       streamline */
				    if (tp->loop_check) {
					int j, index;
					for (j=0; j<stream.npts_ordered; j++) {
					    index = nii_voxel3_index(tp->seed_mask,
								     floor(stream.pts_ordered[3*j+0]),
								     floor(stream.pts_ordered[3*j+1]),
								     floor(stream.pts_ordered[3*j+2]));
					    tp->loop_mask[index] = 0;
					}
				    }
				    
                                    if (stream.npts_ordered > 2) {
                                        track_add_stream_to_buffer(tp, &stream);
					/* write out the seed point for this track. */
					if (tp->save_seeds && fseed_pt != NULL) {
					    fprintf(fseed_pt, "%f,%f,%f,%d\n",
						    stream.seedpt[0], 
						    stream.seedpt[1],
						    stream.seedpt[2],
						    stream.npts_ordered);
					}
                                        tv.n_count++;
                                    }
                                }
                            }
                        }
                    }
                    
                }
            }
        }
        /* I am commenting this out until I can track down a memory bug
	   that only seems to affect really large track files and on
	   64-bit linux (faraday.mbi.ufl.edu).
        fprintf(stderr, "Slice %03.0f of %03d.\n", z, tp->data_dims[2]);
        fflush(stderr);
	*/
    }
    
    if (tp->save_seeds && fseed_pt != NULL) {
	fclose(fseed_pt);
    }
    
    /* See above.
    fprintf(stderr, ".. finished with %d streamlines.\n",tv.n_count);
    */
    
    if (tp->seed_plan == SEED_N_CUBED) free(frac_seeds);
    
    // make sure the buffer gets flushed
    track_flush_buffer(tp);
    
    free(tp->track_buffer);
    nifti_image_free(tp->fa_image);
    
    if (tp->loop_check && tp->loop_mask != NULL)
	free(tp->loop_mask);
    
    // tv.n_count now has the true number of fibers, so we seek to the
    // beginning of the file and rewrite the header
    fseek(tp->output_fp, 0, SEEK_SET);
    iobytes = fwrite(&tv, 1, sizeof(tv), tp->output_fp);
    if (iobytes != sizeof(tv)) {
        fprintf(stderr, "Unable to write to output file.\n");
        perror("Reason");
        exit(1);
    }
    
    fclose(tp->output_fp);
    
    return 0;
}

// convert x,y,z coordinate to single-dimension array index into nii array
inline int nii_voxel3_index(nifti_image *nim,
                            const int x,
                            const int y,
                            const int z)
{
    return (x + y*nim->nx + z*nim->nx*nim->ny);
}

int get_direction_index(track_params *tp, int direction, int index,
                        float *result)
{
    
    if (direction+1 > tp->num_directions) {
        result[0] = 0; result[1] = 0; result[2] = 0;
        return 0;
    }
    
    result[0] = tp->directions[direction][0][index];
    result[1] = tp->directions[direction][1][index];
    result[2] = tp->directions[direction][2][index];
    
    if (result[0] == 0 && result[1] == 0 && result[2] == 0) {
        return 0;
    } else {
        return 1;
    }
    
}

int get_direction_xyz(track_params *tp, int direction, float *pt,
                      float *result)
{
    
    int index;
    int x,y,z;
    
    x = floor(pt[0]);
    y = floor(pt[1]);
    z = floor(pt[2]);
    
    if (direction+1 > tp->num_directions) {
        result[0] = 0; result[1] = 0; result[2] = 0;
        return 0;
    }
    
    index = nii_voxel3_index(tp->direction_imgs[direction], x, y, z);
    result[0] = tp->directions[direction][0][index];
    result[1] = tp->directions[direction][1][index];
    result[2] = tp->directions[direction][2][index];
    if (result[0] == 0 && result[1] == 0 && result[2] == 0) {
        return 0;
    } else {
        return 1;
    }
    
}

int nii_read_directions_from_files(char **file_names, track_params *tp)
{
    
    int i,n;
    nifti_brick_list *nbl;
    nifti_image **nim;
    int blist[3] = {0,1,2};
    int n_directions = tp->num_directions;
    char file[1024];
    
    nbl = malloc(sizeof(nifti_brick_list)*n_directions);
    nim = malloc(sizeof(nifti_image *)*n_directions);
    if (nbl == NULL || nim == NULL) {
        fprintf(stderr, "Unable to allocate memory for nifti data (direction files)\n");
        exit(1);
    }
    
    n=0;
    for (i=0; i<n_directions; i++) {
        sprintf(file,  "%s%c%s", tp->directions_dir, DIRECTORY_SEPARATOR, file_names[i]);
        nim[n] = nifti_image_read_bricks(file, 3, blist, &(nbl[n]));
        if (nim[n] != NULL) { 
            tp->directions[n] = (float **)(nbl[n].bricks);
            n++;
        }
    }
    
    if (n == 0) {
        fprintf(stderr, "*** No valid direction files found.\n");
        exit(1);
    }
    
    tp->num_directions = n;
    tp->direction_imgs = nim;
    tp->direction_brks = nbl;
    tp->voxel_dims[0] = nim[0]->dx;
    tp->voxel_dims[1] = nim[0]->dy;
    tp->voxel_dims[2] = nim[0]->dz;
    
    return 1;
}

nim_compare_result nii_compare_headers (nifti_image *nim1, nifti_image *nim2)
{
    int i,j;
    float tol = 1e-5;
    int result = EQUAL;
    
    if (nim1->ndim != nim2->ndim)
        result |= NUM_DIM_MISMATCH;
    
    for (i=1; i<nim1->ndim; i++) {
        if (nim1->dim[i] != nim2->dim[i]) 
            result |= DIMENSION_MISMATCH;
    }
    
    if (nim1->dx != nim2->dx || nim1->dy != nim2->dy || nim1->dz != nim2->dz)
        result |= GRIDSPACE_MISMATCH;
    
    for (i=0; i<4; i++) {
        for (j=0; j<4; j++) {
            if (abs(nim1->qto_xyz.m[i][j] - nim2->qto_xyz.m[i][j]) > tol) {
                result |= ORIENT_QFORM_MISMATCH;
            } else if (abs(nim1->sto_xyz.m[i][j] - nim2->sto_xyz.m[i][j]) > tol) {
                result |= ORIENT_SFORM_MISMATCH;
            }
        }
    }
    
    return result;
}

float stream_find_best_direction(streamline *s, track_params *tp, float *best)
{
    int i;
    float pdir[3];
    float dp, max_dp;
    //float m1, m2;
    max_dp = 0;
    dp     = 0;
    pdir[0] = 0;
    pdir[1] = 0;
    pdir[2] = 0;
    
    for (i=0; i<tp->num_directions; i++) {
        get_direction_xyz(tp, i, s->currpt, pdir);
        if (pdir[0] == 0 && pdir[1] == 0 && pdir[2] == 0) continue;
        
        dp = (pdir[0]*s->currdir[0] +
              pdir[1]*s->currdir[1] +
              pdir[2]*s->currdir[2]);
        /*
         m1 = sqrt(pdir[0]*pdir[0] +
         pdir[1]*pdir[0] +
         pdir[2]*pdir[0]);
         m2 = sqrt(s->currdir[0]*s->currdir[0] +
         s->currdir[1]*s->currdir[1] +
         s->currdir[2]*s->currdir[2]);
         dp /= m2;
         */	
        if (dp < 0) {
            dp *= -1;
            pdir[0] *= -1; pdir[1] *= -1; pdir[2] *= -1;
        }
        
        if (dp > max_dp) {
            best[0] = pdir[0];
            best[1] = pdir[1];
            best[2] = pdir[2];
            max_dp = dp;
        }	
    }
    
    return max_dp;
    
}

void stream_reorganize_f2b(streamline *s)
{
    
    int i,n;
    
    n=0;
    
    for (i=(s->npts_bw-1); i>0; i--) {
        s->pts_ordered[3*n+0] = s->pts_bw[3*i+0];
        s->pts_ordered[3*n+1] = s->pts_bw[3*i+1];
        s->pts_ordered[3*n+2] = s->pts_bw[3*i+2];
        n++;
    }
	
    for (i=0; i<s->npts_fw; i++) {
        s->pts_ordered[3*n+0] = s->pts_fw[3*i+0];
        s->pts_ordered[3*n+1] = s->pts_fw[3*i+1];
        s->pts_ordered[3*n+2] = s->pts_fw[3*i+2];
        n++;
    }
    
    s->npts_ordered = n;
    
}

int stream_check_loop(track_params *tp, float *pt, int npts)
{
    int index;
    int visited_pt_index;
    int x = floor(pt[0]);
    int y = floor(pt[1]);
    int z = floor(pt[2]);
    
    if (tp->loop_check == 0 || tp->loop_mask == NULL) return 1;
    
    index = nii_voxel3_index(tp->seed_mask, x, y, z);
    visited_pt_index = (tp->loop_mask)[index];

    /*
     * We cause the tracking to stop if the current voxel has been visited before,
     * and is has not been visited in the last 10 steps (due to sub-voxel stepping).
     */
    if (visited_pt_index != 0 && abs(visited_pt_index - npts) > 10) {
	return 0; /* stop the track */
    } else {
	/* place the current step number in the voxel location for comparison
         * later and continue the track.
	 */
	(tp->loop_mask)[index] = npts;
	return 1;
    }
    
}

int stream_check_term_xyz(track_params *tp, float *pt) 
{
    int index;
    int x = floor(pt[0]);
    int y = floor(pt[2]);
    int z = floor(pt[1]);
    
    if (tp->term_mask == NULL) return 1;
    
    index = nii_voxel3_index(tp->term_mask, x, y, z);
    return ((int *) (tp->term_mask->data))[index] > tp->fa_thr ? 1 : 0;
    
}

int stream_check_fa_xyz(track_params *tp, float *pt) 
{
    int index;
    int x = floor(pt[0]);
    int y = floor(pt[1]);
    int z = floor(pt[2]);
    
    if (tp->fa_image == NULL) return 1;
    
    index = nii_voxel3_index(tp->fa_image, x, y, z);
    return ((float *) (tp->fa_image->data))[index] > tp->fa_thr ? 1 : 0;
    
}

int stream_check_fa_index(track_params *tp, int index)
{
    
    if (tp->fa_image == NULL) return 1;
    
    return ((float *) (tp->fa_image->data))[index] >= tp->fa_thr ? 1 : 0;
    
}

int stream_track_fw(track_params *tp, streamline *s)
{
    return stream_track(tp, s, STREAM_DIRECTION_FW);
}

int stream_track_bw(track_params *tp, streamline *s)
{
    return stream_track(tp, s, STREAM_DIRECTION_BW);
}

int stream_track(track_params *tp, streamline *s, int direction) 
{
    int n;
    float pdir[3];
    float *pts;
    int *npts;
    
    if (direction == 0) return 0;
    
    s->currpt[0]  = s->seedpt[0];
    s->currpt[1]  = s->seedpt[1];
    s->currpt[2]  = s->seedpt[2];
    s->currdir[0] = s->seeddir[0];
    s->currdir[1] = s->seeddir[1];
    s->currdir[2] = s->seeddir[2];
    
    if (direction == STREAM_DIRECTION_BW) {
        s->currdir[0] *= -1;
        s->currdir[1] *= -1;
        s->currdir[2] *= -1;
        pts = s->pts_bw;
        npts = &(s->npts_bw);
    } else {
        pts = s->pts_fw;
        npts = &(s->npts_fw);
    }
    
    n = 0;
    pdir[0] = s->currdir[0];
    pdir[1] = s->currdir[1];
    pdir[2] = s->currdir[2];
    
    do {
        
        float dp;
        
        if (stream_check_fa_xyz(tp, s->currpt) == 0) break;
        if (stream_check_term_xyz(tp, s->currpt) == 0) break;	
        if (stream_check_loop(tp, s->currpt, *npts) == 0) break;
	
        s->currpt[0] += (pdir[0]*tp->step_vec[0]);
        s->currpt[1] += (pdir[1]*tp->step_vec[1]);
        s->currpt[2] += (pdir[2]*tp->step_vec[2]);
        
        if ( (s->currpt[0] < 0 || s->currpt[0] >= tp->data_dims[0]) ||
            (s->currpt[1] < 0 || s->currpt[1] >= tp->data_dims[1]) ||
            (s->currpt[2] < 0 || s->currpt[2] >= tp->data_dims[2]) ) break;
	    
        pts[3*(*npts)+0] = (s->currpt[0])*tp->voxel_dims[0];
        pts[3*(*npts)+1] = (s->currpt[1])*tp->voxel_dims[1];
        pts[3*(*npts)+2] = (s->currpt[2])*tp->voxel_dims[2];

	if ((3 * *npts)+3 >= STREAM_BUFSIZE) {
	    /* one more point would kill us. */
	    break;
	}
	
	(*npts)++;
        
        dp = stream_find_best_direction(s, tp, pdir);
        
        if (pdir[0] == 0 && pdir[1] == 0 && pdir[2] == 0) break;
        if (isnan(pdir[0]) || isnan(pdir[1]) || isnan(pdir[2])) break;
        if (dp < tp->cos_angle_thr) break;

        s->currdir[0] = pdir[0];
        s->currdir[1] = pdir[1];
        s->currdir[2] = pdir[2];
		
    } while (1);
    
    
    if (s->prop_direction == STREAM_DIRECTION_FW) {
        s->prop_direction = STREAM_DIRECTION_BW;
        return 1;
    } else {
        s->prop_direction = STREAM_DONE;
        return 0;
    }
    
}

int nii_recast_to_int32 (nifti_image *nim)
{
    int i, nvoxels;
    int *data;
    
    nvoxels = nim->dim[1]*nim->dim[2]*nim->dim[3];
    data = malloc(sizeof(int)*nvoxels);
    if (data == NULL) {
        fprintf(stderr, "Unable to allocate memory for seed mask.\n");
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
                fprintf(stderr, "Nifti datatype not supported: %s\n", 
                        nifti_datatype_string(nim->datatype));
                return 0;
        }
    }
    
    free(nim->data);
    nim->data = (void *)data;
    
    return 1;
}

void track_flush_buffer(track_params *tp)
{
    int iobytes;
    
    iobytes = fwrite(tp->track_buffer, 1, tp->track_buf_used, tp->output_fp);
    if (iobytes != tp->track_buf_used) {
        perror("track_flush_buffer(): Unable to write to output track file");
        exit(1);
    }
    
    tp->track_buf_used = 0;
    
}

void track_add_stream_to_buffer(track_params *tp, streamline *s)
{
    
    int track_bytes = sizeof(float)*3*s->npts_ordered;
    
    if (tp->track_buf_used+4+track_bytes > tp->track_bufsize) {
        track_flush_buffer(tp);
    }
    
    memcpy((tp->track_buffer+tp->track_buf_used),
           &(s->npts_ordered), 4);
    
    memcpy((tp->track_buffer+tp->track_buf_used+4),
           s->pts_ordered, track_bytes);
    
    tp->track_buf_used += track_bytes+4;
    
}

int track_sanity_check(track_params *tp)
{
    
    int testresult, i;
    int finalresult = 1;
    
    //    return 1;
    
    nifti_image *nim0 = tp->direction_imgs[0];
    nifti_image *nimi = NULL;
    
    if (nim0 == NULL) return 0;
    
    /* step 1: check that the nifti parameters of the direction files match */
    /* Note: we rely on these pointers not being NULL. This should not happen
     here, and should have been detected when we read them in. */
    fprintf(stderr, "### Starting sanity check on input files......");
    for (i=1; i<tp->num_directions; i++) {
        
        nim0 = tp->direction_imgs[0];
        nimi = tp->direction_imgs[i];
        
        testresult = nii_compare_headers(nim0, nimi);
        if (testresult & DIMENSION_MISMATCH) {
            fprintf(stderr, 
                    "*** Dimension mismatch found in direction files:\n");
            fprintf(stderr,
                    "    File 00: [ %d, %d, %d ], File %2d: [ %d, %d, %d ]\n",
                    nim0->nx, nim0->ny, nim0->nz, i, nimi->nx, nimi->ny, nimi->nz);
            finalresult = 0;
        }
        
        if (testresult & GRIDSPACE_MISMATCH) {
            fprintf(stderr,
                    "*** Voxel grid spacing mismatch found in direction files:\n");
            fprintf(stderr, 
                    "*** File 00: [ %f, %f, %f ], File %2d: [ %f, %f, %f ]\n",
                    nim0->dx, nim0->dy, nim0->dz, i, nimi->dx, nimi->dy, nimi->dz);
            finalresult = 0;
        }
        
        if (nimi->dim[0] != 4 || (testresult & NUM_DIM_MISMATCH)) {
            fprintf(stderr,
                    "*** Number of dimensions mismatch found in direction files:\n");
            fprintf(stderr, 
                    "*** File 00: %d, File %2d: %d\n",
                    nim0->dim[0], i, nimi->dim[0]);
            finalresult = 0;
        }
        
        if (testresult & (ORIENT_QFORM_MISMATCH|ORIENT_SFORM_MISMATCH)) {
            fprintf(stderr,
                    "### Warning: Orientation mismatch detected in direction files:\n");
            fprintf(stderr,
                    "### Processing will continue, but use caution when analyzing result.\n");
        }
	    
    }
    
    if (finalresult == 0) {
        fprintf(stderr, "*** Please correct discrepancies in direction files.\n");
        fprintf(stderr, "**** retry tracking.\n");
        return finalresult;
    }
    
    /* step 2: check that the orientation and dimensions of seed mask match
     that of the first direction file */
    
    nimi = tp->seed_mask;
    testresult = nii_compare_headers(nim0, nimi);
    
    if (testresult & DIMENSION_MISMATCH) {
        fprintf(stderr, 
                "*** Dimension mismatch found between direction files and seed mask:\n");
        fprintf(stderr,
                "    Direction File: [ %d, %d, %d ], Seed Mask: [ %d, %d, %d ]\n",
                nim0->nx, nim0->ny, nim0->nz, nimi->nx, nimi->ny, nimi->nz);
        finalresult = 0;
    }
    
    if (testresult & GRIDSPACE_MISMATCH) {
        fprintf(stderr,
                "*** Voxel grid spacing mismatch found between direction files");
        fprintf(stderr,
                "and seed mask:\n");
        fprintf(stderr, 
                "*** Direction File: [ %f, %f, %f ], Seed Mask: [ %f, %f, %f ]\n",
                nim0->dx, nim0->dy, nim0->dz, nimi->dx, nimi->dy, nimi->dz);
        finalresult = 0;
    }
    
    if (nimi->dim[0] != 3) {
        fprintf(stderr,
                "*** Incorrect number of dimensions found in seed mask file:\n");
        fprintf(stderr, 
                "*** Should be: 3, Seed Mask: %d\n", nimi->dim[0]);
        finalresult = 0;
    }
    
    if (testresult & (ORIENT_QFORM_MISMATCH|ORIENT_SFORM_MISMATCH)) {
        fprintf(stderr,
                "### Warning: Orientation mismatch detected in direction files:\n");
        fprintf(stderr,
                "### Processing will continue, but use caution when analyzing result.\n");
    }
    
    
    /* step 3: check that the FA mask matches the seed mask in orientation
     and dimensions */
    
    if (tp->fa_image != NULL) {
        nim0 = tp->fa_image;
        testresult = nii_compare_headers(nimi, nim0);
        
        if (testresult & DIMENSION_MISMATCH) {
            fprintf(stderr, 
                    "*** Dimension mismatch found between FA image and seed mask:\n");
            fprintf(stderr,
                    "    FA Image: [ %d, %d, %d ], Seed Mask: [ %d, %d, %d ]\n",
                    nim0->nx, nim0->ny, nim0->nz, nimi->nx, nimi->ny, nimi->nz);
            finalresult = 0;
        }
        
        if (testresult & GRIDSPACE_MISMATCH) {
            fprintf(stderr,
                    "*** Voxel grid spacing mismatch found between FA image and ");
            fprintf(stderr,
                    "and seed mask:\n");
            fprintf(stderr, 
                    "*** FA Image: [ %f, %f, %f ], Seed Mask: [ %f, %f, %f ]\n",
                    nim0->dx, nim0->dy, nim0->dz, nimi->dx, nimi->dy, nimi->dz);
            finalresult = 0;
        }
        
        if (nim0->dim[0] != 3) {
            fprintf(stderr,
                    "*** Incorrect number of dimensions found in FA image file:\n");
            fprintf(stderr, 
                    "*** Should be: 3, FA Image: %d\n", nimi->dim[0]);
            finalresult = 0;
        }
        
        if (testresult & (ORIENT_QFORM_MISMATCH|ORIENT_SFORM_MISMATCH)) {
            fprintf(stderr,
                    "### Warning: Orientation mismatch detected between seed mask");
            fprintf(stderr,
                    "### and FA image. Processing will continue, but use caution");
            fprintf(stderr,
                    "### when analyzing result.\n");
        }
        
    }
    
    /* step 4: check that the termination mask matches the seed mask in
     orientation and dimensions */
    
    if (tp->term_mask != NULL) {
        nim0 = tp->term_mask;
        testresult = nii_compare_headers(nimi, nim0);
        
        if (testresult & DIMENSION_MISMATCH) {
            fprintf(stderr, 
                    "*** Dimension mismatch found between termination and seed mask:\n");
            fprintf(stderr,
                    "    Termination Mask: [ %d, %d, %d ], Seed Mask: [ %d, %d, %d ]\n",
                    nim0->nx, nim0->ny, nim0->nz, nimi->nx, nimi->ny, nimi->nz);
            finalresult = 0;
        }
        
        if (testresult & GRIDSPACE_MISMATCH) {
            fprintf(stderr,
                    "*** Voxel grid spacing mismatch found between termination and ");
            fprintf(stderr,
                    "and seed mask:\n");
            fprintf(stderr, 
                    "*** Termination Mask: [ %f, %f, %f ], Seed Mask: [ %f, %f, %f ]\n",
                    nim0->dx, nim0->dy, nim0->dz, nimi->dx, nimi->dy, nimi->dz);
            finalresult = 0;
        }
        
        if (nim0->dim[0] != 3) {
            fprintf(stderr,
                    "*** Incorrect number of dimensions found in Termination Mask file:\n");
            fprintf(stderr, 
                    "*** Should be: 3, FA Image: %d\n", nimi->dim[0]);
            finalresult = 0;
        }
        
        if (testresult & (ORIENT_QFORM_MISMATCH|ORIENT_SFORM_MISMATCH)) {
            fprintf(stderr,
                    "### Warning: Orientation mismatch detected between seed mask");
            fprintf(stderr,
                    "### and termination mask. Processing will continue, but use ");
            fprintf(stderr,
                    "### caution when analyzing result.\n");
        }
        
    }
    
    return finalresult;
    
}

void track_print_params(track_params *tp)
{
    
    int den = tp->sd_dens;
    
    fprintf(stderr, "Number of direction files........: %d\n", tp->num_directions);
    if (tp->seed_plan == SEED_N_CUBED) {
        fprintf(stderr, "Seeding Density (n^3)............: %d (%d/voxel)\n", den, den*den*den);
    } else {
        fprintf(stderr, "Seeding Density (randomize)......: %d (randomly placed)\n", den);
    }
    fprintf(stderr, "Streamline Step Size (in mm).....: %.3f\n", tp->step_size);
    fprintf(stderr, "Streamline Step Resolution (vxl).: [%.3f, %.3f, %.3f]\n", 
            tp->step_vec[0], tp->step_vec[1], tp->step_vec[2]);
    fprintf(stderr, "Data Voxel Dimensions (in mm)....: [%.3f, %.3f, %.3f]\n",
            tp->seed_mask->dx, tp->seed_mask->dy, tp->seed_mask->dz);
    fprintf(stderr, "Fractional Anisotropy Threshold..: %.3f\n", tp->fa_thr);
    fprintf(stderr, "Turn Angle Threshold (in deg)....: %.3f\n", tp->angle_thr);
    fprintf(stderr, "Minimum Track Length Threshold...: %.3f\n", tp->stream_len_thr);
    fprintf(stderr, "Seed Mask File...................: %s\n", tp->seed_mask_file);
    fprintf(stderr, "FA mask file.....................: %s\n", tp->fa_mask_file);
    fprintf(stderr, "Loop checking (1=on, 0=off)......: %d\n", tp->loop_check);
    fprintf(stderr, "Save seed points (1=on, 0=off)...: %d\n", tp->save_seeds);
    fprintf(stderr, "Track Output File................: %s\n", tp->output_file);
    fprintf(stderr, "Track Buffer Size (in bytes).....: %d\n", tp->track_bufsize);
    
}

void track_print_usage(void)
{
    printf("Usage: track_track [options]\n\n");
    printf("where [options] is given by the following:\n");
    printf("  -o output_file      - a path and filename to save the tractograhy output\n");
    printf("                        Note that if the file exists it will be overwritten.\n");
    printf("                        (default: tracks_out.trk).\n\n");
    printf("  -d directions_dir   - path to directory containing fiber direction\n");
    printf("                        volumes. They should be nifti volumes.\n");
    printf("                        (default: current directory).\n\n");
    printf("  -sdm seed_mask      - a 3D nifti volume that defines the region(s)\n");
    printf("                        to be seeded for tracking\n");
    printf("                        (default: seedmask.nii.gz).\n\n");
    printf("  -sd seed_density    - an integer indicating how many seed points per\n");
    printf("                        voxel. The number of seed points will be\n");
    printf("                        seed_density^3 (default: 1).\n\n");
    printf("  -sp seed plan       - this is the seed plan. It can be either 'n_cubed'\n");
    printf("                        which will arrange the seed points uniformly throughout\n");
    printf("                        each voxel, or 'random' which will arrange the seed\n");
    printf("                        points uniformly *at random* throughout each voxel.\n");
    printf("                        (default: n_cubed).\n\n");
    printf("  -fa                 - the fractional anisotropy threshold. Tracks entering\n");
    printf("                        a voxel with FA below threshold will be stopped.\n");
    printf("                        (default: 0.05).\n\n");
    printf("  -fam                - a 3D nifti volume containing the fractional anisotropy\n");
    printf("                        data.\n");
    printf("                        (default: <none>).\n\n");
    printf("  -ss step_size       - a numeric argument indicating the step size for\n");
    printf("                        streamline propagation in units of voxels.\n");
    printf("                        (default: 0.5).\n\n");
    printf("  -ang threshold      - the turn angle cutoff threshold in degrees. Streamlines \n");
    printf("                        will be stopped if they take a turn greater than this\n");
    printf("                        angle.\n");
    printf("                        (default: 50 degrees).\n\n");
    printf("  -loop               - Turn on loop checking. Streamlines will be stopped if they");
    printf("                        return to a previously visited voxel.\n");
    printf("                        (default: off)\n\n");
    printf("  -so                 - Save the seed points for each track to a file. The file will be");
    printf("                        created in the same location as the output track file, and will.\n");
    printf("                        have the extension '.seeds' (default: off)\n\n");
    
}

