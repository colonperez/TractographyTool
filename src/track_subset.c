#include "trackvis.h"
#include "../../TrackTools_GUI/tt_current_version.h"

int main (int argc, char **argv)
{
    
    char *trackfile_in;
    char *trackfile_out;
    
    FILE *fp_in, *fp_out;
    int skip_increment, tmpsize, nbytes;
    int i, ntracks_out, track_size, ctr;
    
    tv_header tv;
    float *buffer;
    
    if (argc != 4) {
        printf("Usage: track_subset INPUT_FILE SKIP_INCREMENT OUTPUT_FILE\n\n");
        return 1;
    } else {
        trackfile_in   = argv[1];
        skip_increment = atoi(argv[2]);
        trackfile_out  = argv[3];
    }
    
    fprintf(stderr, "track_subset (%d) starting...\n", kTT_CURRENT_VERSION);
    
    ntracks_out = 0;
    
    buffer = malloc(kTV_POINT_BUFSIZE*sizeof(float));
    memset(buffer, 0, kTV_POINT_BUFSIZE*sizeof(float));
    memset(&tv, 0, sizeof(tv));
    
    fp_in = fopen(trackfile_in, "rb");
    
    if (NULL == fp_in) {
        fprintf(stderr, "Unable to open: %s.\n", trackfile_in);
        return 1;
    }
    
    nbytes = fread(&tv, 1, sizeof(tv), fp_in);
    if (nbytes != sizeof(tv) || strcmp(tv.id_string, "TRACK") || tv.hdr_size != 1000) {
        fprintf(stderr, "Bad trackvis file: %s.\n", trackfile_in);
        fclose(fp_in);
        return 1;
    }
    
    if (0 == tv.n_count) {
        tv.n_count = get_track_count(fp_in, &tv);
        if (tv.n_count == -1) {
            fprintf(stderr, "Unable to read %s, cannot determine number of tracks.\n",
                    trackfile_in);
            fclose(fp_in);
            return 1;
        }
    }
    
    printf("Track count in: %d\n", tv.n_count);
    
    fp_out = fopen(trackfile_out, "w+b");
    
    if (NULL == fp_out) {
        fprintf(stderr, "Unable to open %s.\n", trackfile_out);
        return 1;
    }
    
    if (fwrite(&tv, 1, sizeof(tv), fp_out) != sizeof(tv)) {
        fprintf(stderr, "Unable to write trackvis header!\n");
        perror("Reason");
        return 1;
    }
    ctr = 0;
    
    for (i=0; i<tv.n_count; i++) {
        
        fread(&track_size, 1, sizeof(track_size), fp_in);
        
        if (track_size > kTV_POINT_BUFSIZE) {
            fprintf(stderr, "Track %d too large for buffer (track_size = %d), exiting.\n", i,track_size);
            fclose(fp_in);
            return 1;
        }
        
        tmpsize = (track_size * (3 + tv.n_scalars) + tv.n_properties) * sizeof(float);
        fread(buffer, 1, tmpsize, fp_in);
        
        if (0 == i % skip_increment) {
            fwrite(&track_size, 1, sizeof(track_size), fp_out);
            fwrite(buffer, 1, tmpsize, fp_out);
            ntracks_out++;
        }
        
        if (++ctr > ((float) tv.n_count) * 0.02) {
            printf(".");
            fflush(stdout);
            ctr = 0;
        }
        
    }
    
    printf("\nTrack count out: %d\n", ntracks_out);
    tv.n_count = ntracks_out;
    fseek(fp_out, 0, SEEK_SET);
    fwrite(&tv, 1, sizeof(tv), fp_out);
    fclose(fp_in);
    fclose(fp_out);
    
    free(buffer);
    
    return 0;
    
}

