/*
 * sicall.h: Calculate the likelihoods of SNPs and INDELs. Make the calling.
 * Author: Mengyao Zhao
 * Create date: 2011-08-09
 * Contact: zhangmp@bc.edu
 * Last revise: 2013-02-07 
 */

#include <stdlib.h>

/*! @function	Call the SNPs and INDELs based on the trained PHMM parameters.
	@param	ref	reference sequence
	@param	ref_name	reference (chromosome) name
	@param	begin	PHMM start coordinate on the target reference
	@param	filter	SNP and INDEL quality filter 
 */
void likelihood (bamFile fp,
				 bam_index_t* idx,
				double** transition, 
				 double** emission, 
				 char* ref, 
				 int32_t tid, 
				 int32_t window_beg,	// 0_based coordinate
				 int32_t region_beg,	// 0_based coordinate
				 int32_t region_end, 	// 0_based coordinate
				 int32_t size,
				 int32_t filter) {
