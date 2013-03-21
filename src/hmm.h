/*
 * hmm.h: Banded semi-global HMM
 * Author: Mengyao Zhao
 * Create date: 2011-06-13
 * Contact: zhangmp@bc.edu
 * Last revise: 2013-03-21 
 */

#include <stdlib.h>

/*!	@typedef	Structure for read blocks
	@field	count	number of reads
	@field	pos		0-based leftmost coordinate
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
	int32_t* pos;
	int32_t* seq_l;	// read length
	uint8_t* seqs;
} reads;

/*! @function	Initialize the transition matrixes array. */
double** transition_init (const double a, const double b, const double r, const double c, const double d, const int32_t L);

/*! @function	Destroy the transition matrixes array. */
void transition_destroy (double** matrix_array, const int32_t L);

/*! @function	Initialize the emission matrix. */
double** emission_init (char* ref, int32_t size);

/*! @function	Destroy the emission matrix. */
void emission_destroy (double** array, const int32_t L);

/*! @function	Full possibility forward and backward algorithm 
 */
double forward_backward (double** transition, 
						 double** emission,
						 int32_t ref_begin,
						 int32_t window_len,	// region size 
						 uint8_t* read, 
						 int32_t read_len, 
						 double** f, 
						 double** b, 
						 double* s, 
						 int32_t bw);

/*! @function	Baum-Welch algorithm for parameter estimation
	@param	ref_len	reference sequence length 
 */
void baum_welch (double** transition, 
				 double** emission, 
				 char* ref_seq, 
				 int32_t window_begin,	// 0-based coordinate 
				 int32_t window_len,
				 int32_t bw,	// band width 
				 reads* r, 
				 double df);

