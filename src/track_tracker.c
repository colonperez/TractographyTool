
#include "../../TrackTools_GUI/tt_current_version.h"
#include "track_track.h"

void track_initialize_opts(track_params *tp, int argc, char **argv);

int main (int argc, char **argv) {
    
    track_params tp;
    
    memset(&tp, 0, sizeof(tp));
    
    // initialize the track_params using command-line opts
    track_initialize_opts(&tp, argc, argv);
    
    fprintf(stderr, "track_tracker (%d) starting...\n", kTT_CURRENT_VERSION);
    
    if (track_sanity_check(&tp) == 0) {
        fprintf(stderr, "Sanity check failed.\n");
        exit(1);
    } else {
        fprintf(stderr, "Sanity check passed.\n");
    }
    
    track_print_params(&tp);
    
    return track_perform_tracking(&tp);
    
}

void track_initialize_opts(track_params *tp, int argc, char **argv)
{
    int opt;
    float max;
    
    char *direction_files_multi[5] = {
        "V_00_all.nii.gz", 
        "V_01_all.nii.gz",
        "V_02_all.nii.gz",
        "V_03_all.nii.gz",
        "V_04_all.nii.gz"};
    
    tp->fa_thr         = 0.1;
    tp->point_limit    = STREAM_BUFSIZE;
    tp->track_bufsize  = TRACK_BUFSIZE;
    tp->sd_dens        = 1;
    tp->output_file    = "tracks_out.trk";
    tp->fa_thr         = 0.05;
    tp->angle_thr      = 50.0;
    tp->cos_angle_thr  = cos(50.0*M_PI/180.0);
    tp->seed_mask_file = "seedmask.nii.gz";
    tp->directions_dir = ".";
    tp->step_size      = 0.5;
    tp->num_directions = 5;
    tp->seed_plan      = SEED_N_CUBED;
    tp->loop_check     = 0;
    
    for (opt=1; opt<argc; opt++) {
        if (0 == strcmp(argv[opt], "-h")) {
            track_print_usage();
            exit(0);
        }
        if (0 == strcmp(argv[opt], "-d")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -d requires an argument.\n");
                exit(1);
            }
            tp->directions_dir = argv[opt+1];
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-o")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -o requires an argument.\n");
                exit(1);
            }	    
            tp->output_file = argv[opt+1];
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-loop")) {
	    tp->loop_check = 1;
	    continue;
        } else if (0 == strcmp(argv[opt], "-so")) {
	    tp->save_seeds = 1;
	    continue;
        } else if (0 == strcmp(argv[opt], "-fa")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -f requires an argument.\n");
                exit(1);
            }	    
            tp->fa_thr = atof(argv[opt+1]);
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-fam")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -fam requires an argument.\n");
                exit(1);
            }	    
            tp->fa_mask_file = argv[opt+1];
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-sd")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -sd requires an argument.\n");
                exit(1);
            }	    
            tp->sd_dens = atoi(argv[opt+1]);
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-sp")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -sd requires an argument.\n");
                exit(1);
            } else if (0 == strcmp(argv[opt+1], "n_cubed")) {
                tp->seed_plan = SEED_N_CUBED;
            } else if (0 == strcmp(argv[opt+1], "random")) {
                tp->seed_plan = SEED_RANDOM;
            } else {
                fprintf(stderr, "Error: -sp requires either n_cubed or random\n");
                exit(1);
            }
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-ss")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -ss requires an argument.\n");
                exit(1);
            }	    
            tp->step_size = atof(argv[opt+1]);
            opt++;
            continue;	    
        } else if (0 == strcmp(argv[opt], "-tbs")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -tbs requires an argument.\n");
                exit(1);
            }	    
            tp->track_bufsize = atoi(argv[opt+1]);
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-sdm")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -sdm requires an argument.\n");
                exit(1);
            }	    
            tp->seed_mask_file = argv[opt+1];
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-term")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -term requires an argument.\n");
                exit(1);
            }	    
            tp->term_mask_file = argv[opt+1];
            opt++;
            continue;
        } else if (0 == strcmp(argv[opt], "-ang")) {
            if (opt+1 == argc || argv[opt+1][0] == '-') {
                fprintf(stderr, "Error: -ang requires a numeric argument.\n");
                exit(1);
            }	    
            tp->angle_thr     = atof(argv[opt+1]);
            tp->cos_angle_thr = cos(tp->angle_thr*M_PI/180.0);
            opt++;
            continue;
        } else {
            fprintf(stderr, "Ignoring junk on command line: %s\n", argv[opt]);
            fprintf(stderr, "Use -h to see a list of options.\n");
            exit(1);
        }
    }
    
    fprintf(stderr, "Loading direction files...\n");
    nii_read_directions_from_files(direction_files_multi, tp);
    fprintf(stderr, "...done (%d direction files found).\n", tp->num_directions);
    
    fprintf(stderr, "Loading seed mask...");
    tp->seed_mask = nifti_image_read(tp->seed_mask_file, 1);
    if (tp->seed_mask == NULL) { 
        fprintf(stderr, "... unable to load %s\n", tp->seed_mask_file);
        exit(1);
    }
    nii_recast_to_int32(tp->seed_mask);
    fprintf(stderr, " done.\n");
    
    if (tp->term_mask_file != NULL) {
        fprintf(stderr, "Loading termination mask...");
        tp->term_mask = nifti_image_read(tp->term_mask_file, 1);
        if (tp->term_mask == NULL) { 
            fprintf(stderr, "... unable to load %s\n", tp->term_mask_file);
            exit(1);
        }
        nii_recast_to_int32(tp->term_mask);
        fprintf(stderr, " done.\n");
    }
    
    if (tp->fa_thr > 0.0 && tp->fa_mask_file != NULL) {
        fprintf(stderr, "Loading FA mask...");
        tp->fa_image = nifti_image_read(tp->fa_mask_file, 1);
        if (tp->fa_image == NULL) {
            fprintf(stderr, "... unable to load %s\n", tp->fa_mask_file);
            exit(1);
        } else {
            fprintf(stderr, "... done.\n");
        }
    }
    
    fprintf(stderr, "Allocating completed track buffer...");
    tp->track_buffer = malloc(tp->track_bufsize*sizeof(char));
    if (tp->track_buffer == NULL) {
        fprintf(stderr, " unable to allocate buffer.\n");
        exit(1);
    } else {
        fprintf(stderr, " done.\n");
    }
    tp->track_buf_used = 0;
    
    tp->data_dims[0] = tp->seed_mask->nx;
    tp->data_dims[1] = tp->seed_mask->ny;
    tp->data_dims[2] = tp->seed_mask->nz;

    if (tp->loop_check == 1) {
	fprintf(stderr, "Allocating loop-check buffer...");
	tp->loop_mask = malloc(sizeof(int)*(tp->data_dims[0]*tp->data_dims[1]*tp->data_dims[2]));
	if (tp->loop_mask == NULL) {
	    fprintf(stderr, " unable to allocate buffer for loop checking.\n");
	    exit(1);
	} else {
	    fprintf(stderr, " done.\n");
	}
    } else {
	fprintf(stderr, "Loop-check not enabled (use -loop to enable).\n");
    }
    
    max = tp->seed_mask->dx;
    if (tp->seed_mask->dy > max) max = tp->seed_mask->dy;
    if (tp->seed_mask->dz > max) max = tp->seed_mask->dz;
    
    tp->step_vec[0] = tp->seed_mask->dx/max*tp->step_size;
    tp->step_vec[1] = tp->seed_mask->dy/max*tp->step_size;
    tp->step_vec[2] = tp->seed_mask->dz/max*tp->step_size;
    
}

