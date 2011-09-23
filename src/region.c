/*
 * region.c: Get reference and alignments in a region using samtools-0.1.16
 * Author: Mengyao Zhao
 * Create date: 2011-06-05
 * Last revise data: 2011-09-23
 * Contact: zhangmp@bc.edu 
 */

#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

int main (int argc, char * const argv[]) {
	int32_t ret = 0;	//return value
	fprintf (stdout, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

	if (argc == optind) return usage();
	if (argc == optind + 2) { // convert/print the entire file
		fprintf (stderr, "Please give out a genome region with the length shorter or equal to 10M.\n\n");
		fprintf(stderr, "Notes:\n\
\n\
     A region should be presented in one of the following formats:\n\
     `chr1', `chr2:1,000' and `chr3:1000-2,000'. When a region is\n\
     specified, the input alignment file must be an indexed BAM file.\n\
\n");} else {	// retrieve alignments in specified regions
		int32_t i;
		bamFile fp;
		bam_header_t* header;
		bam1_t* bam = bam_init1();
		bam_index_t *idx = 0;
	    if (faidx_t* fai = fai_load(argv[1]) == 0) {
			fprintf(stderr, "Random alignment retrieval requires the reference index file.\n");
			ret = 1;
			goto end_fai;
		}
		if (fp = bam_open(argv[2], "r") == 0) {
			printf(stderr, "Fail to open \"%s\" for reading.\n", argv[2]);
			ret = 1;
			goto end_fp;
		}
		if (header = bam_header_read(fp) == 0) {
			fprintf(stderr, "Fail to read the header from \"%s\".\n", argv[2]);
			ret = 1;
			goto end_header;
		}
		if (idx = bam_index_load(argv[2]) == 0) { // index is unavailable
			fprintf(stderr, "Random alignment retrieval only works for indexed BAM files.\n");
			ret = 1;
			goto end_idx;
		}
		for (i = optind + 2; i < argc; ++i) {
			int32_t tid, beg, end, result, ref_len = 0;
			char* ref_seq = "";
			int32_t n = 70, l = 65536;
			int32_t half_len = 0, count = 0;
			bam_parse_region(header, argv[i], &tid, &beg, &end); // parse a region in the format like `chr2:100-200'
			
			if (tid < 0) { // reference name is not found
				fprintf(stderr, "Region \"%s\" specifies an unknown reference name. Continue anyway.\n", argv[i]);
				ret = 1;
				goto end;
			}
			char* coordinate = ":beg-end";
			char* region = calloc(strlen(header->target_name[tid]) + strlen(coordinate) + 1, sizeof(char));
			strcpy(region, header->target_name[tid]);
			strcat(region, coordinate);
			ref_seq = fai_fetch(fai, region, &ref_len); /* ref_len is a return value */
			if (ref_seq == 0) {
				fprintf(stderr, "Retrieval of reference region \"%s\" failed due to truncated file or corrupt reference index file\n", argv[i]);
				ret = 1;
				goto end;
			}
			free(region);
			
			double** transition = transition_init (0.002, 0.98, 0.00067, 0.02, 0.998, ref_len);
			double** emission = emission_init(ref_seq);

			reads* r = calloc(1, sizeof(reads));
			r->seq_l = calloc(n, sizeof(int32_t));
			r->seqs = calloc(l, sizeof(uint8_t));

			bam_iter_t bam_iter = bam_iter_query(idx, tid, beg, end);
			while (bam_iter_read (fp, bam_iter, bam) >= 0) {
				uint8_t* read_seq = bam1_seq(bam);
				int32_t read_len = bam->core.l_qseq;
				int32_t char_len = read_len/2, j;
				if (count >= n) {
					n = count + 1;
					kroundup32(n);
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
				r->seq_l[count] = read_len;
				count ++;
			}
			r->count = count;

			bam_iter_destroy(bam_iter);

			baum_welch (transition, emission, ref_seq, ref_len, r, 0.01); /* 0-based coordinate */ 

			free(r->seqs);
			free(r->seq_l);
			free(r);

			likelihood (transition, emission, ref_seq, header->target_name[i], 0, 0);	
			emission_destroy(emission, ref_len);
			transition_destroy(transition, ref_len);
			free (ref_seq);
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
	}
	return ret;
}

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   region <in.fasta> <in.bam> [region1 [...]]\n\n");
	fprintf(stderr, "Notes:\n\
\n\
     A region should be presented in one of the following formats:\n\
     `chr1', `chr2:1,000' and `chr3:1000-2,000'. When a region is\n\
     specified, the input alignment file must be an indexed BAM file.\n\
\n");
	return 1;
}
