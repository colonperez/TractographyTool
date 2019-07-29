#include "trackvis.h"
#include "../../TrackTools_GUI/tt_current_version.h"

#define MAXFILES 20
int tv_split_interleaved(int argc, char **argv);

int main (int argc, char **argv)
{
    
    return tv_split_interleaved(argc, argv);
    
}

int tv_split_interleaved (int argc, char **argv) 
{
    
    char *trackfile_in;
    char *trackfile_out_pref;
    char trackfile_out[MAXFILES][256+1];
    
    FILE *fps_out[MAXFILES];
    int ntracks_out[MAXFILES];
    
    int num_files;
    int i, ntracks_in, nbytes;
    
    tv_file *tvf;
    tv_header *tvh;
    tv_track *tv_tk;
    
    if (argc != 4) {
        printf("Usage: track_split INPUT_FILE NUMBER_OF_FILES OUTPUT_FILE_PREFIX\n\n");
        printf("track_split splits the input file into <number_of_files> new files of\n");
        printf("approximately the same size. The maximum number of files is 20.\n");
        printf("Output files will be names <output_file_prefix>_????, where ????\n");
        printf("is a number in a four digit series starting with '0000'.\n\n");
        return 1;
    } else {
        trackfile_in        = argv[1];
        num_files           = atoi(argv[2]);
        if (num_files > MAXFILES) {
            fprintf(stderr, "The maximum number of output files is %d.\n\n", MAXFILES);
            return 1;
        }
        
        trackfile_out_pref  = argv[3];
    }
    
    fprintf(stderr, "track_split (%d) starting...\n", kTT_CURRENT_VERSION);
    
    tvf = tv_open(trackfile_in);
    if (NULL == tvf) {
        fprintf(stderr, "Unable to open: %s.\n", trackfile_in);
        return 1;
    } else {
        tvh = tvf->tv_hdr;
    }
    
    ntracks_in = tvh->n_count;
    printf("Track count in: %d\n", ntracks_in);
    
    for (i=0; i<num_files; i++) {
        sprintf(trackfile_out[i],  "%s_%04d.trk", trackfile_out_pref, i);
        fps_out[i] = fopen(trackfile_out[i], "w+b");
        if (NULL == fps_out[i]) {
            fprintf(stderr, "Unable to open %s.\n", trackfile_out[i]);
            return 1;
        }
        nbytes = fwrite(tvh, 1, sizeof(*tvh), fps_out[i]);
        if (nbytes != sizeof(*tvh)) {
            fprintf(stderr, "Unable to write trackvis header.\n");
            perror("Reason");
            return 1;
        }
    }
    
    tv_tk = tv_new_track();
    
    for (i=0; i<ntracks_in; i++) {
        
        if (tv_read_next_track(tvf, tv_tk) == 0) {
            fprintf(stderr, "Unable to read from track file.\n");
            return 1;
        }
		
        nbytes = fwrite(&(tv_tk->n_points), 1, 4, fps_out[i % num_files]);
        nbytes = fwrite(tv_tk->t_data, 1, tv_tk->n_elements*sizeof(float),
                        fps_out[i % num_files]);
        if (nbytes != tv_tk->n_elements*sizeof(float)) {
            fprintf(stderr, "Unable to write to trackvis file.\n");
            perror("Reason");
            return 1;
        }
        
        ntracks_out[i % num_files]++;
        
    }
    
    for (i=0; i<num_files; i++) {
        tvh->n_count = ntracks_out[i];
        fseek(fps_out[i], 0, SEEK_SET);
        nbytes = fwrite(tvh, 1, sizeof(*tvh), fps_out[i]);
        if (nbytes != sizeof(*tvh)) {
            fprintf(stderr, "Unable to write trackvis header.\n");
            perror("Reason");
            return 1;
        }	
        fclose(fps_out[i]);
    }
    
    tv_close(tvf);
    tv_free_track(tv_tk);
    
    return 0;
    
}

