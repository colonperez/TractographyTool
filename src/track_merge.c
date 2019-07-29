#include "trackvis.h"
#include "../../TrackTools_GUI/tt_current_version.h"

int main (int argc, char **argv)
{
    tv_file *tvf;
    tv_header tvh;
    char *trackfile_out;    
    char *buffer;
    int ntracks_out, f, nbytes;
    FILE *fp_out;
    
    if (argc < 2) {
        printf("Usage: track_merge OUTPUT_FILE INPUT_FILE_1 INPUT_FILE_2 ...\n\n");
        return 1;
    } else {
        trackfile_out = argv[1];
    }
    
    fprintf(stderr, "track_merge (%d) starting...\n", kTT_CURRENT_VERSION);
    
    ntracks_out = 0;
    
    buffer = malloc(sizeof(char)*kTV_POINT_BUFSIZE);
    memset(buffer, 0, kTV_POINT_BUFSIZE);
    memset(&tvh, 0, sizeof(tvh));
    
    fp_out = fopen(trackfile_out, "w+b");
    if (NULL == fp_out) {
        fprintf(stderr, "Unable to open %s.\n", trackfile_out);
        return 1;
    }
    
    nbytes = fwrite(&tvh, 1, sizeof(tvh), fp_out);
    if (nbytes != sizeof(tvh)) {
        fprintf(stderr, "An error occurred trying to write track file header.\n");
        return 1;
    } 
    
    for (f=2; f<argc; f++) {
        
        //fp_in = fopen(argv[f], "rb");
        tvf = tv_open(argv[f]);
        if (NULL == tvf) {
            fprintf(stderr, "Unable to open: %s.\n", argv[f]);
            return 1;
        } else {
            fprintf(stderr, "Reading: %s\n", argv[f]);
        }
        
        fseek(tvf->tv_fp, 1000, SEEK_SET);
        
        while ((nbytes = fread(buffer, 1, kTV_POINT_BUFSIZE, tvf->tv_fp)) != 0) {
            fwrite(buffer, nbytes, 1, fp_out);
        }
        
        ntracks_out += tvf->tv_hdr->n_count;
        
        tv_close(tvf);
        
        
        
    }
    /* FIXME: add error checking to i/o functions */
    tvh.n_count = ntracks_out;
    fseek(fp_out, 0, SEEK_SET);
    fwrite(&tvh, 1, sizeof(tvh), fp_out);
    fclose(fp_out);
    free(buffer);
    
    return 0;
    
}
