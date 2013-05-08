/*
 * sicall.h: Calculate the likelihoods of SNPs and INDELs. Make the calling.
 * Author: Mengyao Zhao
 * Create date: 2011-08-09
 * Contact: zhangmp@bc.edu
 * Last revise: 2013-05-08 
 */

#include <stdlib.h>
#include "bam.h"
#include "khash.h"

#ifndef KHASH
#define KHASH
KHASH_MAP_INIT_INT(insert, char*)
KHASH_MAP_INIT_INT(mnp, char*)
#endif
/*
#ifndef NUM2BASE
#define NUM2BASE
char num2base (int8_t num);
#endif
*/
/*! @function	Call the SNPs and INDELs based on the trained PHMM parameters.
	@param	ref	reference sequence
	@param	ref_name	reference (chromosome) name
	@param	begin	PHMM start coordinate on the target reference
	@param	filter	SNP and INDEL quality filter 
 */
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
				khash_t(mnp) *hm);
