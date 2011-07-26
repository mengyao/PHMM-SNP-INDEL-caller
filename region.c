/*
 * Get reference and alignments in a region using samtools-0.1.16
 * Author: Mengyao Zhao
 * Create date: 2011-06-05
 * Last revise data: 2011-07-25
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
	int len = 0;

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
		char* coordinate = ":0-999";
		char* region = calloc(strlen(header->target_name[i]) + strlen(coordinate) + 1, sizeof(char));
		char* ref_seq;
		int32_t n = 70, l = 65536;
		reads* r = calloc(1, sizeof(reads));
		r->seq_l = calloc(n, sizeof(int32_t));
		r->seqs = calloc(l, sizeof(uint8_t));

		strcpy(region, header->target_name[i]);
		strcat(region, coordinate);
		ref_seq = fai_fetch(fai, region, &len); /* len is a return value */
		free(region);
		
		bam1_t* bam = bam_init1();
		bam_index_t* idx = bam_index_load(argv[2]);
		bam_iter_t bam_iter = bam_iter_query(idx, i, 0, 999);
		int32_t half_len = 0, count = 0;
		while (bam_iter_read (fp, bam_iter, bam) >= 0) {
			uint8_t* read_seq = bam1_seq(bam);
			char* read_name = bam1_qname(bam);
			int32_t read_len = bam->core.l_qseq;
			if (count >= n) {
				kroundup32(n);
				r->seq_l = realloc(r->seq_l, n * sizeof(int32_t));	
			}
			if (half_len * 2 >= l - 8192) {
				kroundup32(l);
				r->seqs = realloc(r->seqs, l * sizeof(uint8_t));
			}		
	
			fprintf (stdout, "read name: %s\n", read_name);
			fprintf (stdout, "read_seq in BW: %s\n", read_seq);
			int32_t char_len = read_len/2, j;
		/*	if (total_len%2) {
				for (j = (total_len + 1)/2; j < (total_len + 1)/2 + char_len; j ++) {
					r->seqs[j - 1] = (r->seqs[j - 1] & 0xf0) | (read_seq[j] >> 4 & 0xf);
					r->seqs[j] = read_seq[j] << 4;
				}
				if (read_len%2) r->seqs[j - 1] = (r->seqs[j - 1] & 0xf0) | (read_seq[j] >> 4 & 0xf);
			} else {*/
			for (j = half_len; j < half_len + char_len; j ++) {
				r->seqs[j] = read_seq[j];
			}
			if (read_len%2) r->seqs[j] = read_seq[j];
		
			half_len += char_len + read_len%2;
			r->seq_l[count] = read_len;
			count ++;
		}
		r->count = count;

		bam_iter_destroy(bam_iter);
		bam_index_destroy(idx);
		bam_destroy1(bam);

		baum_welch (ref_seq, len, r, 0.01); /* 0-based coordinate */ 
		
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
