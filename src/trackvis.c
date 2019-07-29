#include "trackvis.h"

tv_track *tv_new_track(void)
{
    
    tv_track *tvt = malloc(sizeof(tv_track));
    
    if (tvt == NULL) {
        fprintf(stderr, "Unable to create track data structure: no memory\n");
        return tvt;
    }
    
    tvt->t_data = malloc(sizeof(float) * kTV_POINT_BUFSIZE);
    if (tvt->t_data == NULL) {
        fprintf(stderr, "Unable to allocate memory for track point buffer.\n");
        free(tvt);
        tvt = NULL;
        return tvt;
    }
    
    return tvt;
    
}

void tv_free_track(tv_track *tvt)
{
    free(tvt->t_data);
    tvt->t_data = NULL;
    free(tvt);
    tvt = NULL;
}

tv_file *tv_create(const char *filename, tv_header *tvh)
{
    tv_file *tvf = NULL;
    FILE *fp = NULL;
    
    fp = fopen(filename, "w+b");
    if (fp == NULL) {
        fprintf(stderr, "Unable to open or create trackvis file.\n");
        perror("Reason");
        return tvf;
    }
    
    if (tvh == NULL) {
        fprintf(stderr, "tv_create requires a properly allocated tv_header\n");
        return tvf;
    }
    
    tvf = malloc(sizeof(tv_file));
    if (tvf == NULL) {
        fprintf(stderr, "Unable to allocate memory for trackvis header\n");
        return tvf;
    }
    
    /* zero out the count since this is a new file. */
    tvh->n_count = 0;
    
    tvf->tv_hdr = tvh;
    tvf->tv_fp = fp;
    
    return tvf;
}

tv_header *tv_clone_header(tv_header *tvh)
{
    
    tv_header *tvh_clone = NULL;
    
    if (tvh == NULL) return tvh_clone;
    
    tvh_clone = malloc(sizeof(tv_header));
    if (tvh_clone == NULL) {
        fprintf(stderr, "tv_clone_header: unable to allocate memory for header.\n");
        return tvh_clone;
    }
    
    memset(tvh_clone, 0, sizeof(*tvh_clone));
    memcpy(tvh_clone, tvh, sizeof(*tvh));
    
    /* Create default header using parameters from seed volume.
     See: http://www.trackvis.org/docs/?subsect=fileformat */
    /* handle with memcpy()
     tvh_clone->id_string[0] = 'T'; tvh_clone->id_string[1] = 'R';
     tvh_clone->id_string[2] = 'A'; tvh_clone->id_string[3] = 'C';
     tvh_clone->id_string[4] = 'K'; tvh_clone->id_string[5] = '\0';
     tvh_clone->dim[0] = tvh->dim[0];
     tvh_clone->dim[1] = tvh->dim[1];
     tvh_clone->dim[2] = tvh->dim[2];
     tvh_clone->voxel_size[0] = tvh->voxel_size[0];
     tvh_clone->voxel_size[1] = tvh->voxel_size[1];
     tvh_clone->voxel_size[2] = tvh->voxel_size[2];
     tvh_clone->version = tvh->version;
     tvh_clone->hdr_size = tvh->hdr_size;
     tvh_clone->image_orientation_patient[0] = tvh->image_orientation_patient[0];
     tvh_clone->image_orientation_patient[1] = tvh->image_orientation_patient[1];
     tvh_clone->image_orientation_patient[2] = tvh->image_orientation_patient[2];
     tvh_clone->image_orientation_patient[3] = tvh->image_orientation_patient[3];
     tvh_clone->image_orientation_patient[4] = tvh->image_orientation_patient[4];
     tvh_clone->image_orientation_patient[5] = tvh->image_orientation_patient[5];
     tvh_clone->voxel_order[0] = tvh->voxel_order[0];
     tvh_clone->voxel_order[1] = tvh->voxel_order[1]; 
     tvh_clone->voxel_order[2] = tvh->voxel_order[2];
     tvh_clone->n_count = tvh->n_count;
     */
    
    return tvh_clone;
    
}

tv_file *tv_open(const char *filename)
{
    
    tv_file *tvf = NULL;
    tv_header *tvh = NULL;
    int iobytes = 0;
    FILE *fp_in;
    
    fp_in = fopen(filename, "rb");
    
    if (NULL == fp_in) {
        fprintf(stderr, "Unable to open track file.\n");
        perror("Reason");
        return tvf;
    }
    
    tvh = malloc(sizeof(tv_header));
    
    memset(tvh, 0, sizeof(tv_header));
    
    iobytes = fread(tvh, 1, sizeof(tv_header), fp_in);
    
    if (strcmp(tvh->id_string, "TRACK") || tvh->hdr_size != 1000) {
        fprintf(stderr, "Bad trackvis file: %s.\n", filename);
        fclose(fp_in);
        free(tvh);
        return tvf;
    } else if (iobytes != sizeof(tv_header)) {
        fclose(fp_in);
        perror("Reason");
        free(tvh);
        return tvf;
    }
    
    if (0 == tvh->n_count) {
        tvh->n_count = get_track_count(fp_in, tvh);
        if (tvh->n_count == -1) {
            fprintf(stderr, "Unable to read %s, cannot determine number of tracks.\n", 
                    filename);
            fclose(fp_in);
            free(tvh);
            return tvf;
        }
    }
    
    tvf = malloc(sizeof(tv_file));
    tvf->tv_hdr = tvh;
    tvf->tv_fp = fp_in;
    
    return tvf;
    
}

void tv_close(tv_file *tvf)
{
    if (tvf->tv_hdr != NULL) {
        free(tvf->tv_hdr);
        tvf->tv_hdr = NULL;
    }
    if (tvf->tv_fp != NULL) {
        fclose(tvf->tv_fp);
        tvf->tv_fp = NULL;
    }
}


int tv_read_next_track(tv_file *tvf, tv_track *tvt)
{
    
    int iobytes = 0;
    int trksize = 0;
    int tmpsize = 0;
    off_t fpos;
    
    tv_header *tvh = tvf->tv_hdr;
    
    iobytes = fread(&trksize, 1, 4, tvf->tv_fp);
    
    if (iobytes != 4) {
        perror("Error: Unable to read track size from track file");
        return 0;
    } else if (trksize <= 0) {
        fprintf(stderr, "Error: track size unreasonable (track_size = %d).\n",
                trksize); 
    } else if (trksize > kTV_POINT_BUFSIZE) {
	fpos = ftell(tvf->tv_fp);
        fprintf(stderr, "Error: track too large for buffer (track_size = %d, fpod = %ld).\n",
                trksize, fpos);
        return 0;
    }
    
    tmpsize = (trksize * (3 + tvh->n_scalars) + tvh->n_properties) * sizeof(float);
    iobytes = fread(tvt->t_data, 1, tmpsize, tvf->tv_fp);
    
    if (iobytes != tmpsize) {
        perror("Error: Unable to read track data from track file");
        return 0;
    }
    
    tvt->n_points = trksize;
    tvt->n_elements = 3*trksize;
    return 1;
    
}    

int get_track_count(FILE *fp_in, tv_header *tv)
{
    
    int track_size = 0;
    int next_offset = 0;
    int track_count = 0;
    
    fseek(fp_in, tv->hdr_size, SEEK_SET);
    
    while (0 != fread(&track_size, sizeof(track_size), 1, fp_in)) {
        
        if (track_size == 0) {
            fprintf(stderr, "get_track_count: zero track size for %d\n", track_count);
            return -1;
        }
        
        next_offset = (track_size * (3 + tv->n_scalars) + tv->n_properties) * sizeof(float);
        
        fseek(fp_in, next_offset, SEEK_CUR);
        
        track_count++;
        
    }
    
    return track_count;
    
}

float get_track_length(tv_header *tv, float *track, int track_size) {
    
    int j;
    float length, dist;
    float x,y,z;
    float px,py,pz;
    float dx,dy,dz;
    
    px = track[0];
    py = track[3*0+1];
    pz = track[3*0+2];
    
    dist = 0;
    length = 0;
    
    for (j=1; j<track_size; j++) {
        x = track[3*j+0];
        y = track[3*j+1];
        z = track[3*j+2];
        
        dx = px - x;
        dy = py - y;
        dz = pz - z;
        
        dist = dx*dx + dy*dy + dz*dz;
        dist = sqrt(dist);
        length += dist;
        
        px = x;
        py = y;
        pz = z;
    }
    
    return length;
    
}

void init_tv_header_from_nii(tv_header *tv, nifti_image *nii)
{
    int n, i,j,k, ijk[3];
    
    memset(tv, 0, sizeof(*tv));
    
    /* Create default header using parameters from seed volume.
     See: http://www.trackvis.org/docs/?subsect=fileformat */
    tv->id_string[0] = 'T'; tv->id_string[1] = 'R';
    tv->id_string[2] = 'A'; tv->id_string[3] = 'C';
    tv->id_string[4] = 'K'; tv->id_string[5] = '\0';
    tv->dim[0] = nii->nx;
    tv->dim[1] = nii->ny;
    tv->dim[2] = nii->nz;
    tv->voxel_size[0] = nii->dx;
    tv->voxel_size[1] = nii->dy;
    tv->voxel_size[2] = nii->dz;
    tv->version = 1;
    tv->hdr_size = 1000;
    tv->image_orientation_patient[0] = 1;
    tv->image_orientation_patient[1] = 0;
    tv->image_orientation_patient[2] = 0;
    tv->image_orientation_patient[3] = 0;
    tv->image_orientation_patient[4] = 1;
    tv->image_orientation_patient[5] = 0;
    tv->voxel_order[0] = 'L';
    tv->voxel_order[1] = 'P'; 
    tv->voxel_order[2] = 'S';
    tv->n_count = 0;
    
    ijk[0] = 0; 
    ijk[1] = 0; 
    ijk[2] = 0;
    nifti_mat44_to_orientation( nii->qto_xyz , &i,&j,&k );
    ijk[0] = i;
    ijk[1] = j;
    ijk[2] = k;
    
    for (n=0; n<3; n++) {
        switch (ijk[n]) {
            case NIFTI_L2R: tv->voxel_order[n] = 'R'; break;
            case NIFTI_R2L: tv->voxel_order[n] = 'L'; break;
            case NIFTI_P2A: tv->voxel_order[n] = 'A'; break;
            case NIFTI_A2P: tv->voxel_order[n] = 'P'; break;
            case NIFTI_I2S: tv->voxel_order[n] = 'S'; break;
            case NIFTI_S2I: tv->voxel_order[n] = 'I'; break;
        }
    }
    
}
