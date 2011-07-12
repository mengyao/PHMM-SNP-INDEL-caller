/*
 * Get reference and alignments in a region using samtools-0.1.16
 * Author: Mengyao Zhao
 * Create date: 2011-06-05
 * Last revise data: 2011-07-12
 * Contact: zhangmp@bc.edu 
 */

#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bam.h"
#include "faidx.h"
#include "hmm.h"

int main (int argc, char * const argv[]) {
	/*
	  Declarations for reference file *.fa
	 */
	int32_t i;
	faidx_t* fai = fai_load(argv[1]);
	int len = 0;

	/*
	  Declarations for bam file *.bam
	 */	
	bamFile fp = bam_open(argv[2], "r");
	bam_header_t* header;
	header = bam_header_read(fp);
	
	/*
	 Declarations for bam index *.bai
	 */
	for (i = 0; i < header->n_targets; i ++) {
		char* coordinate = ":0-999";
		char* region = calloc(strlen(header->target_name[i]) + strlen(coordinate) + 1, sizeof(char));
		char* ref_seq;

		strcpy(region, header->target_name[i]);
		strcat(region, coordinate);
		ref_seq = fai_fetch(fai, region, &len); /* len is a return value */
		free(region);
		
		baum_welch (ref_seq, len, fp, argv[2], i, 0, 999); /* 0-based coordinate */ 
		free(ref_seq);
	}
	bam_header_destroy(header);
	bam_close(fp);
	fai_destroy(fai);
	return 0;
}
