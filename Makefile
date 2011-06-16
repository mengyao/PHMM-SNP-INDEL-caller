CC=			gcc
CFLAGS=		-g -Wall -O2 #-m64 #-arch ppc
DFLAGS=     -D_USE_KNETFILE
LOBJS=		knetfile.o sam_header.o bgzf.o kstring.o bam_aux.o bam.o bam_import.o sam.o bam_index.o \
            bam_pileup.o bam_lpileup.o bam_md.o razf.o faidx.o \
            kprobaln.o
PROG=		region
LIBPATH=

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

.PHONY:all clean cleanlocal
region:$(LOBJS) region.o hmm.o
		$(CC) $(CFLAGS) -o $@ $(LOBJS) region.o hmm.o $(LIBPATH) -lm -lz
sam_header.o:sam_header.h khash.h
phase.o:bam.h khash.h
bamtk.o:bam.h
bam_md.o:bam.h faidx.h
cleanlocal:
		rm -fr *.o a.out *.exe $(PROG) *~ 

clean:cleanlocal-recur
