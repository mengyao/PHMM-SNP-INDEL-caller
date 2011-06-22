/*
 * Get reference and alignments in a region using samtools-0.1.16
 * Author: Mengyao Zhao
 * Create date: 2011-06-05
 * Last revise data: 2011-06-13
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
	bam1_t* b = bam_init1();
	header = bam_header_read(fp);
	
	/*
	 Declarations for bam index *.bai
	 */
	bam_index_t* idx = bam_index_load(argv[2]);
	bam_iter_t bam_iter;

	for (i = 0; i < header->n_targets; i ++) {
		char* coordinate = ":0-999";
		char* region = calloc(strlen(header->target_name[i]) + strlen(coordinate) + 1, sizeof(char));
		char* ref_seq;
		float** matrix_array;
		int32_t m;
		int32_t n;

		strcpy(region, header->target_name[i]);
		strcat(region, coordinate);
		ref_seq = fai_fetch(fai, region, &len); /* len is a return value */
		free(region);

		/* fprintf (stdout, "reference sequence: %s\n", ref_seq); */
		matrix_array = transition_init (0.3, 0.5, 0.2, 0.5, 0.5, len);
		for (m = 0; m < len + 1; m ++) {
			for (n = 0; n < 11; n ++) {
				fprintf (stdout, "%f\t", matrix_array[m][n]);
			}
			fprintf (stdout, "\n");
		}		


		bam_iter = bam_iter_query(idx, i, 0, 999);
 		while (bam_iter_read (fp, bam_iter, b) >= 0) {
			char* read_seq = (char*)bam1_seq(b);
			char* read_name = bam1_qname(b);
		/*	fprintf (stdout, "read name: %s\n", read_name); */
		}
		bam_iter_destroy(bam_iter);
		transition_destroy(matrix_array, len);
		free(ref_seq);
	}
	bam_destroy1(b);
	bam_index_destroy(idx);
	bam_header_destroy(header);
	bam_close(fp);
	fai_destroy(fai);
	return 0;
}
