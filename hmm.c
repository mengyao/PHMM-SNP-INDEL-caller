/*
 * Banded semi-global HMM
 * Author: Mengyao Zhao
 * Date: 2011-06-05
 * Contact: zhangmp@bc.edu 
 */

#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include "bam.h"

int main (int argc, char * const argv[]) {
	bamFile fp = bam_open(argv[1], "r");
	bam_header_t* header;
	bam1_t* b;
	header = bam_header_read(fp);
	b = bam_init1();
	while (bam_read1(fp, b) >= 0) {
		char *s = bam_format1(header, b);
		puts(s);
		free(s);
	}
	bam_destroy1(b);
	bam_header_destroy(header);
	bam_close(fp);
	return 0;
}
