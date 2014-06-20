CC=			clang #gcc
CFLAGS=		-g -Wall #-D VERBOSE_DEBUG #-O3 #-m64  
DFLAGS=     -D_USE_KNETFILE
LOBJS=		knetfile.o sam_header.o bgzf.o kstring.o bam_aux.o bam.o bam_import.o sam.o bam_index.o \
            bam_pileup.o bam_lpileup.o bam_md.o razf.o faidx.o \
            kprobaln.o 
PROG=		phv
LIBPATH=

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

.PHONY:all clean cleanlocal
phv:$(LOBJS) phv.o region.o hmm.o sicall.o viterbi.o repeat.o
		$(CC) $(CFLAGS) -o $@ $(LOBJS) phv.o region.o hmm.o sicall.o viterbi.o $(LIBPATH) -lm -lz
sam_header.o:sam_header.h khash.h
phase.o:bam.h khash.h
bamtk.o:bam.h
bam_md.o:bam.h faidx.h
hmm.o:hmm.h
sicall.o:sicall.h
viterbi.o:viterbi.h
region.o:region.h
cleanlocal:
		rm -fr *.o a.out *.exe $(PROG) *~ 

clean:cleanlocal
