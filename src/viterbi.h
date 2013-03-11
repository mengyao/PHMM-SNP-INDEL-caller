/*
 * viterbi.h: Viterbi algorithm used to generate the insert sequence, genotype MNPs and generate a better set of alignments
 * Author: Mengyao Zhao
 * Create date: 2011-12-12
 * Contact: zhangmp@bc.edu
 * Last revise: 2013-03-04 
 */

#include <stdlib.h>
/*
typedef struct {
	char* seq;
	int32_t pos;
} p_insters;

typedef struct {
	p_inserts* ins;
	uint32_t* cigar;	
	int32_t cigarLen;	
} p_viterbi;
*/
char* viterbi (double** transition, double** emission, int32_t ref_len, uint8_t* read, int32_t read_len, int32_t bpoint);
