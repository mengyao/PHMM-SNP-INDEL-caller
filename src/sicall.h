/*
 * sicall.h: Calculate the likelihoods of SNPs and INDELs. Make the calling.
 * Author: Mengyao Zhao
 * Create date: 2011-08-09
 * Contact: zhangmp@bc.edu
 * Last revise: 2013-07-17 
 */

#include <stdlib.h>
#include "bam.h"
#include "khash.h"
#include "kstring.h"

#ifndef KHASH
#define KHASH
KHASH_MAP_INIT_INT(insert, kstring_t)
KHASH_MAP_INIT_INT(mnp, kstring_t)
KHASH_MAP_INIT_INT(delet, kstring_t)
#endif

void likelihood (bam_header_t* header,
				double** transition, 
				 double** emission, 
				 char* ref,
				 uint16_t* depth, 
				 int32_t tid, 
				 int32_t window_beg,	// 0_based coordinate
				 int32_t region_beg,	// 0_based coordinate
				 int32_t region_end, 	// 0_based coordinate
				 int32_t size,
				 int32_t filter,
				khash_t(insert) *hi,
				khash_t(mnp) *hm,
				khash_t(delet) *hd);
