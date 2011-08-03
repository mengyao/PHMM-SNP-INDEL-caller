/*
 * Get reference and alignments in a region using samtools-0.1.16
 * Author: Mengyao Zhao
 * Create date: 2011-06-05
 * Last revise data: 2011-07-28
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
	int32_t i/*, j*/;
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
		char* coordinate = ":0-999";
		char* region = calloc(strlen(header->target_name[i]) + strlen(coordinate) + 1, sizeof(char));
		char* ref_seq;
		int32_t n = 70, l = 65536/*, total_hl = 0*/;
		reads* r = calloc(1, sizeof(reads));
		r->seq_l = calloc(n, sizeof(int32_t));
		r->seqs = calloc(l, sizeof(uint8_t));

		strcpy(region, header->target_name[i]);
		strcat(region, coordinate);
		ref_seq = fai_fetch(fai, region, &ref_len); /* len is a return value */
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
				n = count + 1;
				kroundup32(n);
fprintf (stderr, "n: %d\tcount: %d\n", n, count);
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

/* test begin*/
/*		float** transition = transition_init (0.3, 0.5, 0.2, 0.5, 0.5, ref_len);
		float** emission = emission_init(ref_seq);

		for (j = 0; j < r->count; j ++) {
			fprintf (stderr, "j: %d\n", j);
			uint8_t* read_seq = &r->seqs[total_hl];
			total_hl += r->seq_l[j]/2 + r->seq_l[j]%2;
			int32_t read_len = r->seq_l[j];
			fb* f = (fb*)calloc(1, sizeof(fb));
			fb* b = (fb*)calloc(1, sizeof(fb));
			f->match = (double**)calloc(read_len, sizeof(double*));
			f->insertion = (double**)calloc(read_len, sizeof(double*));
			f->deletion = (double**)calloc(read_len, sizeof(double*));
			b->match = (double**)calloc(read_len, sizeof(double*));
			b->insertion = (double**)calloc(read_len, sizeof(double*));
			b->deletion = (double**)calloc(read_len, sizeof(double*));
			for (i = 0; i < read_len; i ++) {
				f->match[i] = (double*)calloc(ref_len + 1, sizeof(double));
				f->insertion[i] = (double*)calloc(ref_len + 1, sizeof(double));
				f->deletion[i] = (double*)calloc(ref_len + 1, sizeof(double));
				b->match[i] = (double*)calloc(ref_len + 1, sizeof(double));
				b->insertion[i] = (double*)calloc(ref_len + 1, sizeof(double));
				b->deletion[i] = (double*)calloc(ref_len + 1, sizeof(double));
			}	
			double* s = (double*)calloc(read_len + 1, sizeof(double));
			forward_backward (transition, emission, ref_seq, read_seq, read_len, f, b, s);
			
			free(s);
			for (i = 0; i < read_len; i ++) {
				free(f->match[i]);
				free(f->insertion[i]);
				free(f->deletion[i]);
				free(b->match[i]);
				free(b->insertion[i]);
				free(b->deletion[i]);
			}	
			free(f);
			free(b);
		}
		emission_destroy(emission, ref_len);
		transition_destroy(transition, ref_len);
*/
/* test end */

		baum_welch (ref_seq, ref_len, r, 0.001); /* 0-based coordinate */ 
		
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
