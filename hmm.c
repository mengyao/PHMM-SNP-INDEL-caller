/*
 * Banded semi-global HMM
 * Author: Mengyao Zhao
 * Create date: 2011-06-13
 * Contact: zhangmp@bc.edu 
 */

#include "hmm.h"

const float emission[15][15] = {
  /* A     C     G     T     K     M     R     Y     S     B     V     H     D     N     X */
	{1,    0,    0,    0,    0,    0.5,  0.5,  0,    0,    0,    0.33, 0.33, 0.33, 0.25, 0},	/*A*/
	{0,    1,    0,    0,    0,    0.5,  0,    0.5,  0.5,  0.33, 0.33, 0.33, 0,    0.25, 0},	/*C*/
	{0,    0,    1,    0,    0.5,  0,    0.5,  0,    0.5,  0.33, 0.33, 0,    0.33, 0.25, 0},	/*G*/
	{0,    0,    0,    1,    0.5,  0,    0,    0.5,  0,    0.33, 0,    0.33, 0.33, 0.25, 0},	/*T*/
	{0,    0,    0.5,  0.5,  1,    0,    0.25, 0.25, 0.25, 0.7,  0.17, 0.17, 0.7,  0.5,  0},	/*K*/
	{0.5,  0.5,  0,    0,    0,    1,    0.25, 0.25, 0.25, 0.17, 0.7,  0.7,  0.17, 0.5,  0},	/*M*/
	{0.5,  0,    0.5,  0,    0.25, 0.25, 1,    0,    0.25, 0.17, 0.7,  0.17, 0.7,  0.5,  0},	/*R*/
	{0,    0.5,  0,    0.5,  0.25, 0.25, 0,    1,    0.25, 0.7,  0.17, 0.7,  0.17, 0.5,  0},	/*Y*/
	{0,    0.5,  0.5,  0,    0.25, 0.25, 0.25, 0.25, 1,    0.7,  0.7,  0.17, 0.17, 0.5,  0},	/*S*/
	{0,    0.33, 0.33, 0.33, 0.7,  0.17, 0.17, 0.7,  0.7,  1,    0.49, 0.49, 0.49, 0.75, 0},	/*B*/
	{0.33, 0.33, 0.33, 0,    0.17, 0.7,  0.7,  0.17, 0.7,  0.49, 1,    0.49, 0.49, 0.75, 0},	/*V*/
	{0.33, 0.33, 0,    0.33, 0.17, 0.7,  0.17, 0.7,  0.17, 0.49, 0.49, 1,    0.49, 0.75, 0},	/*H*/
	{0.33, 0,    0.33, 0.33, 0.7,  0.17, 0.7,  0.17, 0.17, 0.49, 0.49, 0.49, 1,    0.75, 0},	/*D*/
	{0.25, 0.25, 0.25, 0.25, 0.5,  0.5,  0.5,  0.5,  0.5,  0.75, 0.75, 0.75, 0.75, 1,    0},	/*N*/
	{0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    1},	/*X*/
};

float** transition_init (const float a, const float b, const float r, const float c, const float d, const int32_t L)
{
	float** matrix_array = calloc (L + 1, sizeof(float*));
	int32_t i;
	for (i = 0; i < L + 1; i ++) {
		matrix_array[i] = calloc (11, sizeof(float));
	}

	/* k = 0: inseart before the reference */	
	for (i = 0; i < 11; i ++) {
		matrix_array[0][i] = 0;
	}
	matrix_array[0][4] = (1 - c)*r;	/* I_k -> M_k+1 */
	matrix_array[0][5] = c*r;	/* I_k -> I_k */
	matrix_array[0][6] = r;	/* I_k -> E */
	matrix_array[0][10] = (1 - d)/(L + 1);	/* S -> I_k */

	/* k = 1 ... L - 1 */
	for (i = 1; i < L; i ++) {
		/* Sum of the following 4 lines equals to 1. */
		matrix_array[i][0] = (1 - 2*a)*(1 - r);	/* M_k -> M_k+1 */
		matrix_array[i][1] = a*(1 - r);	/* M_k -> I_k */
		matrix_array[i][2] = a*(1 - r);	/* M_k -> D_k+1 */
		matrix_array[i][3] = r;	/* M_k -> E */

		/* Sum of the following 3 lines equals to 1. */
		matrix_array[i][4] = (1 - c)*r;	/* I_k -> M_k+1 */
		matrix_array[i][5] = c*r;	/* I_k -> I_k */
		matrix_array[i][6] = r;	/* I_k -> E */

		/* Sum of the following 2 lines equals to 1. */
		matrix_array[i][7] = b;	/* D_k -> M_k+1 */
		matrix_array[i][8] = 1 - b;		/* D_k -> D_k+1 */
		
		/* Sum of the following 2 lines equals to 1. */
		matrix_array[i][9] = d/L;	/* S -> M_k+1 */
		matrix_array[i][10] = (1 - d)/(L + 1);	/* S -> I_k */
	}
	matrix_array[1][7] = 0;	/* State D_1 doesn't exist. */	

	/* k = L */
	for (i = 0; i < 11; i ++) {
		matrix_array[L][i] = 0;
	}
	matrix_array[L][5] = 1 - r;	/* I_k -> I_k */
	matrix_array[L][6] = r;	/* I_k -> E */
	matrix_array[L][10] = (1 - d)/(L + 1);	/* S -> I_k */
	
	return matrix_array;
} 

void transition_destroy (float** matrix_array, const int32_t L)
{
	int32_t i;
	for (i = 0; i < L + 1; i ++) {
		free(matrix_array[i]);
	}
	free(matrix_array);
}

/*
float forward (float** transition, float** emission, char* ref, char* read)
{
	int32_t i;	 iter of read 
	int32_t k;	 iter of reference 
	int32_t ref_len = strlen(ref);
	int32_t read_len = 2*strlen(read);
	
	float f_0_s = 1;
	for (k = 1; k <= ref_len; k ++) {
		
	}
}
*/
