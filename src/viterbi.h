/*
 * viterbi.h: Viterbi algorithm used to generate the insert sequence, genotype MNPs and generate a better set of alignments
 * Author: Mengyao Zhao
 * Create date: 2011-12-12
 * Contact: zhangmp@bc.edu
 * Last revise: 2014-08-12 
 */

#include <stdlib.h>
#include "hmm.h"
#include "khash.h"
#include "kstring.h"

#ifndef KHASH
#define KHASH
KHASH_MAP_INIT_INT(insert, kstring_t)
KHASH_MAP_INIT_INT(mnp, kstring_t)
KHASH_MAP_INIT_INT(delet, kstring_t)
#endif

//extern const int8_t nt_table[128];
extern char num2base[16];

// Generate the hash of insertion, mnp and deletion
void hash_imd (double** transition, 
				double** emission, 				 
				char* ref_seq, 
				 int32_t window_begin,	// 0-based coordinate 
				 int32_t window_len, 
				 int32_t bw,
				 reads* r,
				khash_t(insert) *hi,	// key: 1-based relative position in window; value: insert_str1,insert_str2... (insert_str1 == insert_str2 is possible)
				khash_t(mnp) *hm,
				khash_t(delet) *hd);
