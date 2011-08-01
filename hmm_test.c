/*
 * This program is used to test hmm.c.
 * Author: Mengyao Zhao
 * Create date: 2011-06-28
 * Contact: zhangmp@bc.edu 
 */

#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kseq.h"
#include "hmm.h"
KSEQ_INIT(gzFile, gzread)

int main (int argc, char * const argv[]) {
	// validate argument count
	if (argc < 2) {
		fprintf (stderr, "USAGE: %s <ref.fa/fq> <reads.fa/fq>\n", argv[0]);
		exit (1);
	}
	
	gzFile ref_fp;
	kseq_t *ref_seq;
	int l;
	int m;
	ref_fp = gzopen(argv[1], "r");
	ref_seq = kseq_init(ref_fp);
	while ((l = kseq_read(ref_seq)) >= 0) {
		printf("ref_seq: %s\n", ref_seq->seq.s);
		gzFile read_fp = gzopen(argv[2], "r");
		kseq_t*	read_seq = kseq_init(read_fp);
		int32_t refLen = strlen(ref_seq->seq.s); 
		float** matrix_array = transition_init (0.3, 0.5, 0.2, 0.5, 0.5, refLen);
		float** emission = emission_init(ref_seq->seq.s);
		int32_t i;
		int32_t j;
		fprintf(stdout, "transition matrix:\n");
		for (i = 0; i <= refLen; i ++) {
			for (j = 0; j < 11; j ++) {
				fprintf(stdout, "%f\t", matrix_array[i][j]);
			}
			fprintf(stdout, "\n");
		}
		fprintf(stdout, "emission matrix:\n");
		for (i = 0; i < refLen; i ++) {
			for (j = 0; j < 16; j ++) {
				fprintf(stdout, "%f\t", emission[i][j]);
			}
			fprintf(stdout, "\n");
		}
		while ((m = kseq_read(read_seq)) >= 0) {
			printf("read_seq: %s\n", read_seq->seq.s); 
			
	/*		float forward_score = forward (matrix_array, emission, ref_seq->seq.s, read_seq->seq.s);
			float backward_score = backward (matrix_array, emission, ref_seq->seq.s, read_seq->seq.s);
			fprintf(stdout, "forward: %f\tbackward: %f\n", forward_score, backward_score);
			if (forward_score == 0) {

			} else {}
	*/
			int32_t ref_len = strlen(ref_seq->seq.s);
			int32_t read_len = strlen(read_seq->seq.s);
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
			
			forward_backward(matrix_array, emission, ref_seq->seq.s, (uint8_t*)read_seq->seq.s, read_len, f, b, s);	
	
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
		emission_destroy(emission, refLen);
		transition_destroy(matrix_array, refLen);
		kseq_destroy(read_seq);
		gzclose(read_fp);
	}
	kseq_destroy(ref_seq);
	gzclose(ref_fp);
	return 0;
}
