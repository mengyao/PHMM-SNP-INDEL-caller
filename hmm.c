/*
 * Banded semi-global HMM
 * Author: Mengyao Zhao
 * Date: 2011-06-05
 * Contact: zhangmp@bc.edu 
 */

#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bam.h"
#include "faidx.h"

int main (int argc, char * const argv[]) {
	/*
	  Declarations for reference file *.fa
	 */
	int32_t i;
	faidx_t* fai = fai_load(argv[1]);
	char* ref_seq;
	int len = 0;

	/*
	  Declarations for bam file *.bam
	 */	
	bamFile fp = bam_open(argv[2], "r");
	bam_header_t* header;
	bam1_t* b = bam_init1();
	header = bam_header_read(fp);
	
	/*
	 Declarations for bam index *.bai
	 */
	bam_index_t* idx = bam_index_load(argv[2]);
	bam_iter_t bam_iter;

	fprintf (stdout, "Reference sequences:\n");
	for (i = 0; i < header->n_targets; i ++) {
		char* coordinate = ":0-999";
		char* region = strcat(header->target_name[i], coordinate);
		ref_seq = fai_fetch(fai, region, &len); /* len is a return value */
		fprintf (stdout, "%s\n", ref_seq);
		bam_iter = bam_iter_query(idx, i, 0, 999);
 		while (bam_iter_read (fp, bam_iter, b) >= 0) {
			char* read_seq = (char*)bam1_seq(b);
			char* read_name = bam1_qname(b);
			fprintf (stdout, "read name: %s\n", read_name);
		}
		bam_iter_destroy(bam_iter);
	}
	bam_destroy1(b);
	bam_index_destroy(idx);
	bam_header_destroy(header);
	bam_close(fp);
	return 0;
}
