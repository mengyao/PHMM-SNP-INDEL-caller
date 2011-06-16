/*
 * Banded semi-global HMM
 * Author: Mengyao Zhao
 * Create date: 2011-06-13
 * Contact: zhangmp@bc.edu 
 */

#include "hmm.h"

float** transition_init (const float a, const float b, const float r, const float c, const int32_t L, float* last)
{
	float** matrix_array = calloc (L - 1, sizeof(float*));
	int32_t i;
	for (i = 0; i < L - 1; i ++) {
		matrix_array[i] = calloc (18, sizeof(float));
	}	

	/*	     	 M_k,                  I_k,           D_k,         S,      E       */
	/*	M_k-1	 [0] (1 - 2*A)*(1 -r), [1] a*(1 - r), [2] a*(1 - r), [3] 0,  [4] r,  */
	/*	I_k-1	 [5] (1 - b)*(1 - r),  [6] b*(1 - r), [7] 0,         [8] 0,  [9] r,  */
	/*	D_k-1	 [10] 1 - b,           [11] 0,        [12] b,        [13] 0, [14] 0, */
	/*  M_k                              */
	/*  I_k*/
	matrix_array[0] = (1 - 2*a)*(1 - r);
	matrix_array[1] = matrix_array[2] = a*(1 - r);
	matrix_array[5] = (1 - b)*(1 - r);
	matrix_array[6] = b*(1 - r);
	matrix_array[10] = 1 - b;
	matrix_array[4] = matrix_array[9] = r;
	matrix_array[12] = b;
	matrix_array[3] = matrix_array[7] = matrix_array[8] = matrix_array[11] = matrix_array[13] = matrix_array[14] = 0;

	/*   	M_k, 			I_k       */
	/*	S	[15] (1 - a)/L, [16] a/L, */
    /*  E	[17] 0,         [18] 0,   */

	float** matrix = calloc(5, sizeof(float*));
	for (i = 0; i < 5; i ++) {
		matrix[i] = calloc(5, sizeof(float));
	}
	
	matrix[0][0] = (1 - 2*a)*(1 - r);
	matrix[1][0] = (1 - b)*(1 - r);
	matrix[1][1] = b*(1 - r);
	matrix[2][0] = 1 - b;
	matrix[2][2] = b;
	matrix[3][0] = (1 - c)/L;
	matrix[3][1] = c/L;
	matrix[0][1] = matrix[0][2] = a*(1 - r);
	matrix[0][4] = matrix[1][4] = r;
	matrix[0][3] = matrix[1][2] = matrix[1][3] = matrix[2][1] = matrix[2][3] = matrix[2][4] = matrix[3][2] =\
	matrix[3][3] = matrix[3][4] = matrix[4][0] = matrix[4][1] = matrix[4][2] = matrix[4][3] = matrix[4][4] = 0;
	return matrix;
} 

void transition_destroy (float** matrix_array)
{
	int32_t i;
	for (i = 0; i < 5; i ++) {
		free(matrix[i]);
	}
	free(matrix);
}

float forward (float** transition, float** emission, char* ref, char* read)
{
	int32_t i;	/* iter of read */
	int32_t k;	/* iter of reference */
	int32_t ref_len = strlen(ref);
	int32_t read_len = 2*strlen(read);
	
	float f_0_s = 1;
	for (k = 1; k <= ref_len; k ++) {
		m
	}
}

