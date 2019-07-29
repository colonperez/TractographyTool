track_merge.o: track_merge.c trackvis.h nifti1_io.h nifti1.h znzlib.h \
  ../../TrackTools_GUI/tt_current_version.h
track_intersect.o: track_intersect.c \
  ../../TrackTools_GUI/tt_current_version.h trackvis.h nifti1_io.h \
  nifti1.h znzlib.h
track_split.o: track_split.c trackvis.h nifti1_io.h nifti1.h znzlib.h \
  ../../TrackTools_GUI/tt_current_version.h
track_network.o: track_network.c \
  ../../TrackTools_GUI/tt_current_version.h track_network.h trackvis.h \
  nifti1_io.h nifti1.h znzlib.h
track_track.o: track_track.c track_track.h trackvis.h nifti1_io.h \
  nifti1.h znzlib.h knuthrand.h
track_subset.o: track_subset.c trackvis.h nifti1_io.h nifti1.h znzlib.h \
  ../../TrackTools_GUI/tt_current_version.h
track_tracker.o: track_tracker.c \
  ../../TrackTools_GUI/tt_current_version.h track_track.h trackvis.h \
  nifti1_io.h nifti1.h znzlib.h
trackvis.o: trackvis.c trackvis.h nifti1_io.h nifti1.h znzlib.h
mow_recon.o: mow_recon.c ../../TrackTools_GUI/tt_current_version.h \
  mow_recon.h nifti1_io.h nifti1.h znzlib.h track_track.h trackvis.h \
  matrices.h lssq.h nnls.h
nnls.o: nnls.c nnls.h
lssq.o: lssq.c lssq.h matrices.h
matrices.o: matrices.c matrices.h
nifti1_io.o: nifti1_io.c nifti1_io.h nifti1.h znzlib.h
znzlib.o: znzlib.c znzlib.h
knuthrand.o: knuthrand.c knuthrand.h
track_density_mask.o: track_density_mask.c trackvis.h nifti1_io.h \
  nifti1.h znzlib.h ../../TrackTools_GUI/tt_current_version.h
track_niftimath.o: track_niftimath.c \
  ../../TrackTools_GUI/tt_current_version.h track_niftimath.h \
  nifti1_io.h nifti1.h znzlib.h
