/*
 * Banded semi-global HMM
 * Author: Mengyao Zhao
 * Create date: 2011-06-13
 * Contact: zhangmp@bc.edu 
 */

#include <stdlib.h>

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
	{0.33, 0.33  0,    0.33, 0.17, 0.7,  0.17, 0.7,  0.17, 0.49, 0.49, 1,    0.49, 0.75, 0},	/*H*/
	{0.33, 0,    0.33, 0.33, 0.7,  0.17, 0.7,  0.17, 0.17, 0.49, 0.49, 0.49, 1,    0.75, 0},	/*D*/
	{0.25, 0.25, 0.25, 0.25, 0.5,  0.5,  0.5,  0.5,  0.5,  0.75, 0.75, 0.75, 0.75, 1,    0},	/*N*/
	{0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    1},	/*X*/
};

typedef struct {
	float* matrix1; /* M_k, I_k, D_k -> M_k+1, I_k, D_k+1, S, E */
	float* matrix2; /* S, E -> M_k, I_k */
} transition;

/*! @function	Initialize the transition matrixes array. 
  @param last	pointer to the transition matrix of the last states: I_L, S, E -> I_L, M_L
 */
float** transition_init (const float a, const float b, const float r, const float c, const int32_t L, float* last);

/* Destroy the transition matrixes array. */
void transition_destroy (float** matrix_array);
