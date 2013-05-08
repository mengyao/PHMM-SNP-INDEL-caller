/*
 * viterbi.h: Viterbi algorithm used to generate the insert sequence, genotype MNPs and generate a better set of alignments
 * Author: Mengyao Zhao
 * Create date: 2011-12-12
 * Contact: zhangmp@bc.edu
 * Last revise: 2013-05-08 
 */

#include <stdlib.h>
#include "hmm.h"
#include "khash.h"

#ifndef KHASH
#define KHASH
KHASH_MAP_INIT_INT(insert, char*)
KHASH_MAP_INIT_INT(mnp, char*)
#endif

//#ifndef NUM2BASE
//#define NUM2BASE
char num2base (int8_t num);
//#endif

void hash_insert_mnp (double** transition, 
				double** emission, 				 
				char* ref_seq, 
				 int32_t window_begin,	// 0-based coordinate 
				 int32_t window_len, 
				 int32_t bw, 
				 reads* r,
				khash_t(insert) *hi,
				khash_t(mnp) *hm);
