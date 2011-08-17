/*
 * region.c: Get reference and alignments in a region using samtools-0.1.16
 * Author: Mengyao Zhao
 * Create date: 2011-06-05
 * Last revise data: 2011-08-16
 * Contact: zhangmp@bc.edu 
 */

#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bam.h"
#include "faidx.h"
#include "hmm.h"

/*! @function
  @abstract  Round an integer to the next closest power-2 integer.
  @param  x  integer to be rounded (in place)
  @discussion x will be modified.
 */
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

int main (int argc, char * const argv[]) {
	/*
	  Declarations for reference file *.fa
	 */
	int32_t i;
	faidx_t* fai = fai_load(argv[1]);
	int ref_len = 0;

	/*
	  Declarations for bam file *.bam
	 */	
	bamFile fp = bam_open(argv[2], "r");
	bam_header_t* header;
	header = bam_header_read(fp);
	
	/*
	 Declarations for bam index *.bai
	 */
	for (i = 0; i < header->n_targets; i ++) {
		int32_t j, k;
		char* coordinate = ":0-999";
		char* region = calloc(strlen(header->target_name[i]) + strlen(coordinate) + 1, sizeof(char));
		char* ref_seq;
		int32_t n = 70, l = 65536;
		reads* r = calloc(1, sizeof(reads));
		r->seq_l = calloc(n, sizeof(int32_t));
		r->seqs = calloc(l, sizeof(uint8_t));

		strcpy(region, header->target_name[i]);
		strcat(region, coordinate);
		ref_seq = fai_fetch(fai, region, &ref_len); /* ref_len is a return value */
		free(region);
		
		double** transition = transition_init (0.002, 0.98, 0.00067, 0.02, 0.998, ref_len);
		double** emission = emission_init(ref_seq);
		bam1_t* bam = bam_init1();
		bam_index_t* idx = bam_index_load(argv[2]);
		bam_iter_t bam_iter = bam_iter_query(idx, i, 0, 999);
		int32_t half_len = 0, count = 0;
		while (bam_iter_read (fp, bam_iter, bam) >= 0) {
			uint8_t* read_seq = bam1_seq(bam);
			char* read_name = bam1_qname(bam);
			int32_t read_len = bam->core.l_qseq;

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
	
			fprintf (stdout, "read name: %s\n", read_name);
			int32_t char_len = read_len/2, j;
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
		bam_index_destroy(idx);
		bam_destroy1(bam);

		baum_welch (transition, emission, ref_seq, ref_len, r, 0.01); /* 0-based coordinate */ 
		for (k = 0; k <= ref_len; k ++) {
			for (j = 0; j < 11; j ++) {
				fprintf (stderr, "t[%d][%d]: %g\t", k, j, transition[k][j]);
			}
			fprintf (stderr, "\n");
		}

		for (k = 0; k <= ref_len; k ++) {
			fprintf (stderr, "em[%d][1]: %g\tem[%d][2]: %g\tem[%d][4]: %g\tem[%d][8]: %g\tem[%d][15]: %g\n", k, emission[k][1], k, emission[k][2], k, emission[k][4], k, emission[k][8], k, emission[k][15]);
			fprintf (stderr, "em[%d][0]: %g\tem[%d][3]: %g\tem[%d][5]: %g\tem[%d][9]: %g\tem[%d][14]: %g\n", k, emission[k][0], k, emission[k][3], k, emission[k][5], k, emission[k][9], k, emission[k][14]);
		}
			
		emission_destroy(emission, ref_len);
		transition_destroy(transition, ref_len);
		free(r->seqs);
		free(r->seq_l);
		free(r);
		free(ref_seq);
	}
	bam_header_destroy(header);
	bam_close(fp);
	fai_destroy(fai);
	return 0;
}
