/*
 * viterbi.h: Viterbi algorithm used to generate the insert sequence and a better set of alignments
 * Author: Mengyao Zhao
 * Create date: 2011-12-12
 * Contact: zhangmp@bc.edu
 * Last revise: 2012-05-18 
 */

#include <stdlib.h>

typedef struct {
	char* seq;
	int32_t pos;
} p_insters;

typedef struct {
	p_inserts* ins;
	uint32_t* cigar;	
	int32_t cigarLen;	
} p_viterbi;

p_viterbi* viterbi (double** transition, double** emission, int32_t ref_len, uint8_t* read, int32_t read_len);
