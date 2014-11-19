/*
 * phv.c: a SNP and INDEL caller based on PHMM
 * Author: Mengyao Zhao
 * Create date: 2014-06-20
 * Last revise date: 2014-06-20
 * Contact: zhangmp@bc.edu 
 */

#include <time.h>
#include "bam.h"
#include "faidx.h"
#include "region.h"

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   region [options] <reference.fasta> <alignment.bam> [region1 [...]]\n\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "\t-s N\tN is the largest detectable INDEL length. [default: 100]\n");
	fprintf(stderr, "Notes:\n\
\n\
     A region should be presented in one of the following formats:\n\
     `chr1', `chr2:1,000' and `chr3:1,000-2,000'. When a region is\n\
     specified, the input alignment file must be an indexed BAM file.\n\
\n");
	return 1;
}

int main (int argc, char * const argv[]) {
	int32_t ret = 0;	// return value
	float cpu_time;
	clock_t start, end;
	start = clock();

	int32_t l,i;
	int32_t size = 101;	// default largest detectable INDEL length
	bamFile fp;
	bam_header_t* header;
	bam1_t* bam = bam_init1();
	bam_index_t* idx = 0;
	faidx_t* fai;

	fprintf (stdout, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

	// Parse command line.
	while ((l = getopt(argc, argv, "s:")) >= 0) {
		switch (l) {
			case 's': size = atoi(optarg); break;
		}
	}

	if (optind + 2 > argc) return usage();

	if ((fai = fai_load(argv[optind])) == 0) {
		fprintf(stderr, "Random alignment retrieval requires the reference index file.\n");
		ret = 1;
		goto end_fai;
	}
	if ((fp = bam_open(argv[optind + 1], "r")) == 0) {
		fprintf(stderr, "Fail to open \"%s\" for reading.\n", argv[3]);
		ret = 1;
		goto end_fp;
	}
	if ((header = bam_header_read(fp)) == 0) {
		fprintf(stderr, "Fail to read the header from \"%s\".\n", argv[4]);
		ret = 1;
		goto end_header;
	}
	if ((idx = bam_index_load(argv[optind + 1])) == 0) { // index is unavailable
		fprintf(stderr, "Random alignment retrieval only works for indexed BAM files.\n");
		ret = 1;
		goto end_idx;
	}

	if (argc == (optind + 2)) slide_window_whole(fai, fp, header, bam, idx, size); 	// No region is given by the command line.
	else {	// Regions are given by the command line.
		i = optind + 2;
		while(i < argc) {
			int32_t tid, region_begin, region_end;
			bam_parse_region(header, argv[i], &tid, &region_begin, &region_end); // parse a region in the format like `chr2:100-200'
			if (tid < 0) { // reference name is not found
				fprintf(stderr, "region \"%s\" specifies an unknown reference name.\n", argv[i]);
				return 0;
			}
			slide_window_region(fai, fp, bam, idx, header, tid, region_begin, region_end, size);
			++i;
		}
	}

	bam_index_destroy(idx);
end_idx:
	bam_header_destroy(header);
end_header:
	bam_close(fp);
end_fp:
	fai_destroy(fai);
end_fai:
	bam_destroy1(bam);

	end = clock();
	cpu_time = ((float) (end - start)) / CLOCKS_PER_SEC;
	fprintf(stdout, "\n\nCPU time: %f seconds\n", cpu_time);

	return ret;
}
