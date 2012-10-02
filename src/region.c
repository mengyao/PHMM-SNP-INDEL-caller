/*
 * region.c: Get reference and alignments in a region using samtools-0.1.16
 * Author: Mengyao Zhao
 * Create date: 2011-06-05
 * Last revise data: 2012-10-02
 * Contact: zhangmp@bc.edu 
 */

#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "bam.h"
#include "faidx.h"
#include "hmm.h"
#include "sicall.h"

/*! @function
  @abstract  Round an integer to the next closest power-2 integer.
  @param  x  integer to be rounded (in place)
  @discussion x will be modified.
 */
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

typedef struct {
	char* ref_seq;
	int32_t ref_len;
	double** transition;
	double** emission;
} profile;

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   region [options] <in.fasta> <in.bam> [region1 [...]]\n\n");
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

profile* train (faidx_t* fai, 
		   		int32_t tid,	// reference ID
		   		char* ref_name, 
		   		bamFile fp,
		   		bam1_t* bam, 
		   		bam_index_t* idx, 
		   		int32_t beg, 
		   		int32_t size) {	// maximal detectable INDEL size

		int32_t n = 70, l = 65536, half_len = 0, count = 0, end = beg + 999;
		profile* hmm = (profile*)malloc(sizeof(profile));;

		reads* r = calloc(1, sizeof(reads));
		r->pos = calloc(n, sizeof(int32_t));
		r->seq_l = calloc(n, sizeof(int32_t));
		r->seqs = calloc(l, sizeof(uint8_t));
	
		hmm->ref_seq = faidx_fetch_seq(fai, ref_name, beg, end, &hmm->ref_len);	/* region_name and ref_len are return values */

		if (hmm->ref_seq == 0 || hmm->ref_len < 1) {
		fprintf(stderr, "tid: %d\tbeg: %d\tend: %d\n", tid, beg, end);
			fprintf(stderr, "Retrieval of reference region \"%s:%d-%d\" failed due to truncated file or corrupt reference index file\n", ref_name, beg, end);
			hmm = NULL;
		}

		bam_iter_t bam_iter = bam_iter_query(idx, tid, beg, end);	// 1st read mapping position is beg
		while (bam_iter_read (fp, bam_iter, bam) >= 0) {
			uint8_t* read_seq = bam1_seq(bam);
			int32_t read_len = bam->core.l_qseq;
			int32_t char_len = read_len/2, j;
			if (count >= n) {
				n = count + 1;
				kroundup32(n);
				r->pos = realloc(r->pos, n * sizeof(int32_t));	
				r->seq_l = realloc(r->seq_l, n * sizeof(int32_t));	
			}
			if (half_len * 2 >= l - 8192) {
				kroundup32(l);
				l = half_len * 2 + 1;
				r->seqs = realloc(r->seqs, l * sizeof(uint8_t));
			}		
	
			for (j = half_len; j < half_len + char_len; j ++) {
				r->seqs[j] = read_seq[j - half_len];
			}
			if (read_len%2) r->seqs[j] = read_seq[j - half_len];
		
			half_len += char_len + read_len%2;
			r->pos[count] = bam->core.pos;
			r->seq_l[count] = read_len;
			count ++;
		}

		if (count == 0) {
			fprintf (stderr, "There is no read in the given region \"%s:%d-%d\".\n", ref_name, beg, end);
			hmm = NULL;
		} else {
			hmm->transition = transition_init (0.002, 0.98, 0.00067, 0.02, 0.998, hmm->ref_len + size);
			hmm->emission = emission_init(hmm->ref_seq, size);
			r->count = count;
			baum_welch (hmm->transition, hmm->emission, hmm->ref_seq, beg, hmm->ref_len + size, size, r, 0.01); /* 0-based coordinate */ 
		}

		bam_iter_destroy(bam_iter);
		free(r->seqs);
		free(r->seq_l);
		free(r->pos);
		free(r);
		
		return hmm;
}

int main (int argc, char * const argv[]) {
	int32_t ret = 0;	// return value

	float cpu_time;
	clock_t start, end;
	start = clock();

	int32_t l,i;
	int32_t size = 100;	// default largest detectable INDEL length
	bamFile fp;
	bam_header_t* header;
	bam1_t* bam = bam_init1();
	bam_index_t *idx = 0;
	faidx_t* fai;
	profile* hmm;

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
		fprintf(stderr, "Fail to open \"%s\" for reading.\n", argv[2]);
		ret = 1;
		goto end_fp;
	}
	if ((header = bam_header_read(fp)) == 0) {
		fprintf(stderr, "Fail to read the header from \"%s\".\n", argv[2]);
		ret = 1;
		goto end_header;
	}
	if ((idx = bam_index_load(argv[optind + 1])) == 0) { // index is unavailable
		fprintf(stderr, "Random alignment retrieval only works for indexed BAM files.\n");
		ret = 1;
		goto end_idx;
	}
	for (i = optind + 2; i < argc; ++i) {
		int32_t tid, beg, end;
		bam_parse_region(header, argv[i], &tid, &beg, &end); // parse a region in the format like `chr2:100-200'
		/*
		header: pointer to the header structure
		tid: the returned chromosome ID
		beg: the returned start coordinate
		end: the returned end coordinate
		return: 0 on suceess; -1 on failure
			*/

		if (tid < 0) { // reference name is not found
			fprintf(stderr, "Region \"%s\" specifies an unknown reference name.\n", argv[i]);
			ret = 1;
			goto end;
		}

//FIXME: start to make the sliding window
	
		hmm = train(fai, tid, header->target_name[tid], fp, bam, idx, beg, size);
		if (! hmm) fprintf(stderr, "The HMM profile trainning is failed.\n");
		else {
			likelihood (hmm->transition, hmm->emission, hmm->ref_seq, header->target_name[tid], beg, beg, end, 0);	
			
			free(hmm->ref_seq);
			transition_destroy(hmm->transition, hmm->ref_len + size);
			emission_destroy(hmm->emission, hmm->ref_len + size);
		}
		free(hmm);

//FIXME: end of sliding window loop

	}
end:
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

