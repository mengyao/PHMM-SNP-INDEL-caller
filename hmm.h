/*
 * Banded semi-global HMM
 * Author: Mengyao Zhao
 * Create date: 2011-06-13
 * Contact: zhangmp@bc.edu 
 */

#include <stdlib.h>

/*! @function	Initialize the transition matrixes array. */
float** transition_init (const float a, const float b, const float r, const float c, const float d, const int32_t L);

/* Destroy the transition matrixes array. */
void transition_destroy (float** matrix_array, const int32_t L);
