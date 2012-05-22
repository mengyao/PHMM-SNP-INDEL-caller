/*
 * viterbi.c: Generate the insert sequence and the best sequence alignment.
 * Author: Mengyao Zhao
 * Create date: 2012-05-17
 * Contact: zhangmp@bc.edu
 * Last revise: 2012-05-21 
 */

#include "bam.h"
#include "hmm.h"

p_viterbi* viterbi (double** transition, double** emission, int32_t ref_len, uint8_t* read, int32_t read_len) {
	int32_t i;	 /* iter of read */ 
	int32_t k;	 /* iter of reference */
	double temp_m, temp_i, temp_d;

	fb* v = (fb*)calloc(1, sizeof(fb));
	v->match = (double**)calloc(read_len, sizeof(double*));
	v->insertion = (double**)calloc(read_len, sizeof(double*));
	v->deletion = (double**)calloc(read_len, sizeof(double*));
	for (i = 0; i < read_len; i ++) {
		v->match[i] = (double*)calloc(ref_len + 1, sizeof(double));
		v->insertion[i] = (double*)calloc(ref_len + 1, sizeof(double));
		v->deletion[i] = (double*)calloc(ref_len + 1, sizeof(double));
	}

	for (i = 0; i < read_len; ++i) {
		// k = 0
		temp_m = v->match[][] 
		v->insertion[i][0] = emission[0][bam1_seqi(read, i)];
	}
	
	for (i = 0; i < read_len; i ++) {
		free(v->match[i]);
		free(v->insertion[i]);
		free(v->deletion[i]);
	}
	free(v->match);
	free(v->insertion);
	free(v->deletion);	
	free(v);
}



