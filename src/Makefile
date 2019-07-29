
DEBUG=0
CC=gcc
AR=ar
RANLIB=ranlib
CP=cp
RM=rm

#INCFLAGS=-I../../nifticlib-1.1.0/include -I../../zlib-1.2.5
INCFLAGS=-I. -I../../zlib-1.2.5

#LIBS=-L../../nifticlib-1.1.0/lib -lniftiio -lznz -lm -L../../zlib-1.2.5 -lz
LIBS=-L../lib -lniftiio -lznz -lm -L../../zlib-1.2.5 -lz

SRCS=track_merge.c track_intersect.c track_split.c \
	 track_network.c track_track.c track_subset.c \
	 track_tracker.c trackvis.c mow_recon.c nnls.c \
         lssq.c matrices.c nifti1_io.c znzlib.c knuthrand.c \
	 track_density_mask.c track_niftimath.c

LDFLAGS=

ifeq ($(DEBUG),1)
CFLAGS=-O0 -DHAVE_ZLIB -Wall -g 
else

ifeq ($(ARCH),WIN32)
CFLAGS=-O3 -DHAVE_ZLIB -Wall
else
CFLAGS=-O3 -D_LARGEFILE_SOURCE=1 -D_FILE_OFFSET_BITS=64 -DHAVE_ZLIB -Wall -msse2 -ftree-vectorize -ffast-math
endif

endif

DEPFLAGS=-MM

.SUFFIXES: .c .o

.c.o:
	$(CC) -c $(CFLAGS) $(INCFLAGS) $<

all: mow_recon track_merge track_split track_subset track_intersect \
	track_network track_tracker libtracker track_density_mask   \
	track_niftimath

clean:
	$(RM) -f *.o core *.a
	-$(RM) -f depend.mk

realclean: clean
	rm -f *.o core
	rm -rf ../bin/* ../lib/*.a

depend:	
	-$(RM) -f depend.mk
	$(MAKE) depend.mk

depend.mk:
	$(CC) $(DEPFLAGS) $(INCFLAGS) $(SRCS) >> depend.mk

mow_recon: libniftiio.a mow_recon.o nnls.o lssq.o matrices.o libtracker.a
	$(CC) $(LDFLAGS) -o ../bin/mow_recon mow_recon.o nnls.o lssq.o matrices.o $(LIBS) -L../lib -ltracker

track_density_mask: libniftiio.a track_density_mask.o trackvis.o
	$(CC) $(LDFLAGS) -o ../bin/track_density_mask track_density_mask.o trackvis.o $(LIBS)

track_niftimath: libniftiio.a track_niftimath.o 
	$(CC) $(LDFLAGS) -o ../bin/track_niftimath track_niftimath.o $(LIBS)

track_merge: libniftiio.a track_merge.o trackvis.o
	$(CC) $(LDFLAGS) -o ../bin/track_merge track_merge.o trackvis.o $(LIBS)

track_split: libniftiio.a track_split.o trackvis.o
	$(CC) $(LDFLAGS) -o ../bin/track_split track_split.o trackvis.o $(LIBS)

track_subset: libniftiio.a track_subset.o trackvis.o
	$(CC) $(LDFLAGS) -o ../bin/track_subset track_subset.o trackvis.o $(LIBS)

track_intersect: libniftiio.a track_intersect.o trackvis.o
	$(CC) $(LDFLAGS) -o ../bin/track_intersect track_intersect.o trackvis.o $(LIBS) 

track_network: libniftiio.a track_network.o trackvis.o
	$(CC) $(LDFLAGS) -o ../bin/track_network track_network.o trackvis.o $(LIBS) 

track_tracker: libniftiio.a knuthrand.o track_tracker.o track_track.o trackvis.o
	$(CC) $(LDFLAGS) -o ../bin/track_tracker knuthrand.o track_track.o track_tracker.o trackvis.o $(LIBS)

libtracker: libtracker.a

libtracker.a: track_track.o trackvis.o knuthrand.o
	$(AR) -cru ../lib/libtracker.a knuthrand.o track_track.o trackvis.o
	$(RANLIB) ../lib/libtracker.a
	$(CP) *.h ../include/

libniftiio.a: znzlib.a nifti1_io.o
	$(AR) -cru ../lib/libniftiio.a nifti1_io.o 
	$(RANLIB) ../lib/libniftiio.a
	$(CP) nifti1_io.h nifti1.h ../include/

znzlib.a: znzlib.o
	$(AR) -cru ../lib/libznz.a znzlib.o
	$(RANLIB) ../lib/libznz.a
	$(CP) znzlib.h ../include/

include depend.mk
