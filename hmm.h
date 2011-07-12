/*
 * hmm.h: Banded semi-global HMM
 * Author: Mengyao Zhao
 * Create date: 2011-06-13
 * Contact: zhangmp@bc.edu
 * Last revise: 2011-07-12 
 */

#include <stdlib.h>

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
void forward_backward (float** transition, float** emission, char* ref, char* read);

/*! @function	Baum-Welch algorithm for parameter estimation
	@param	len	reference sequence length 
 */
void baum_welch (char* ref_seq, int len, bamFile fp, char* bai, int32_t target_num, int32_t begin, int32_t end); /* 0-based coordinate */ 

