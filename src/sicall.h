/*
 * sicall.h: Calculate the likelihoods of SNPs and INDELs. Make the calling.
 * Author: Mengyao Zhao
 * Create date: 2011-08-09
 * Contact: zhangmp@bc.edu
 * Last revise: 2011-08-24 
 */

#include <stdlib.h>

void likelihood (double** transition, double** emission, char* ref, char* target_name, int32_t begin, int32_t filter);

