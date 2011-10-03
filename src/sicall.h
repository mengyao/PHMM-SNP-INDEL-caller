/*
 * sicall.h: Calculate the likelihoods of SNPs and INDELs. Make the calling.
 * Author: Mengyao Zhao
 * Create date: 2011-08-09
 * Contact: zhangmp@bc.edu
 * Last revise: 2011-10-02 
 */

#include <stdlib.h>

/*! @function	Call the SNPs and INDELs based on the trained PHMM parameters.
	@param	ref	reference sequence
	@param	target_name	reference (chromosome) ID
	@param	begin	PHMM start coordinate on the target reference
	@param	filter	SNP and INDEL quality filter 
 */
void likelihood (double** transition, double** emission, char* ref, char* target_name, int32_t begin, int32_t filter);

