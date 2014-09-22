/*
 * pileup.h: Evaluate the regions with candidate variations from the pileup data using samtools-0.1.18
 * Author: Mengyao Zhao
 * Create date: 2014-07-23
 * Last revise date: 2014-07-25
 * Contact: zhangmp@bc.edu 
 */

#include "bam.h"

int pileup_check (bamFile fp, 
			 		   bam_index_t* idx,
					char* ref_seq, 
			           int32_t tid, 
			  	       int32_t beg,
				   	int32_t end);
