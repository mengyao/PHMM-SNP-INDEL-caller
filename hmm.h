/*
 * hmm.h: Banded semi-global HMM
 * Author: Mengyao Zhao
 * Create date: 2011-06-13
 * Contact: zhangmp@bc.edu
 * Last revise: 2011-06-24 
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

/*! @function	Full possibility forward algorithm 
	@return	Sum of the forward paths probabilities
 */
void forward_backward (float** transition, float** emission, char* ref, char* read);

/*! @function	Full possibility backward algorithm 
	@return	Sum of the backward paths probabilities
 */
/*float backward(float** transition, float** emission, char* ref, char* read);
*/
