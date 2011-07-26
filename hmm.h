/*
 * hmm.h: Banded semi-global HMM
 * Author: Mengyao Zhao
 * Create date: 2011-06-13
 * Contact: zhangmp@bc.edu
 * Last revise: 2011-07-25 
 */

#include <stdlib.h>

/*!	@typedef	Structure for forward and backward matrices
 */
typedef struct {
	double** match;
	double** insertion;
	double** deletion;
	double final;
} fb;

/*!	@typedef	Structure for read blocks
	@field	count	number of reads
	@field	seq_hl	read length
	@field	seqs	read sequences, reads always be stored from the very left side of one byte memory.
					For the reads with length of odd numbers, the lower 4 bits of the last byte memory
					for this read is not used.

	@discussion Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G,
  				8 for T and 15 for N. Two bases are packed in one byte with the base
  				at the higher 4 bits having smaller coordinate on the read. It is
  				recommended to use bam1_seqi() macro to get the base.
 */
typedef struct {
	int32_t count;
	int32_t* seq_l;
	uint8_t* seqs;
} reads;

/*! @function	Initialize the transition matrixes array. */
float** transition_init (const float a, const float b, const float r, const float c, const float d, const int32_t L);

/*! @function	Destroy the transition matrixes array. */
void transition_destroy (float** matrix_array, const int32_t L);

/*! @function	Initialize the emission matrix. */
float** emission_init (char* ref);

/*! @function	Destroy the emission matrix. */
void emission_destroy (float** array, const int32_t L);

/*! @function	Full possibility forward and backward algorithm 
 */
double forward_backward (float** transition, float** emission, char* ref, uint8_t* read, int32_t read_len, fb* f, fb* b);

/*! @function	Baum-Welch algorithm for parameter estimation
	@param	len	reference sequence length 
 */
void baum_welch (char* ref_seq, int ref_len, reads* r, float df); /* 0-based coordinate */ 

