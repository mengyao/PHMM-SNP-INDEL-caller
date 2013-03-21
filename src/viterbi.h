/*
 * viterbi.h: Viterbi algorithm used to generate the insert sequence, genotype MNPs and generate a better set of alignments
 * Author: Mengyao Zhao
 * Create date: 2011-12-12
 * Contact: zhangmp@bc.edu
 * Last revise: 2013-03-21 
 */

#include <stdlib.h>
void hash_insert_mnp (double** transition, 
				double** emission, 				 
				char* ref_seq, 
				 int32_t window_begin,	// 0-based coordinate 
				 int32_t window_len, 
				 int32_t bw, 
				 reads* r,
				khash_t(insert) *hi,
				khash_t(mnp) *hm);
