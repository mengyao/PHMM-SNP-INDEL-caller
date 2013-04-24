CC=			gcc #clang
CFLAGS=		-g -Wall #-O2 #-D VERBOSE_DEBUG #-m64 #-arch ppc 
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
region:$(LOBJS) region.o hmm.o sicall.o viterbi.o
		$(CC) $(CFLAGS) -o $@ $(LOBJS) region.o hmm.o sicall.o viterbi.o $(LIBPATH) -lm -lz
sam_header.o:sam_header.h khash.h
phase.o:bam.h khash.h
bamtk.o:bam.h
bam_md.o:bam.h faidx.h
hmm.o:hmm.h
sicall.o:sicall.h
viterbi.o:viterbi.h
cleanlocal:
		rm -fr *.o a.out *.exe $(PROG) *~ 

clean:cleanlocal-recur
