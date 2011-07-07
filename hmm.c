/*
 * hmm.c: Banded semi-global HMM
 * Author: Mengyao Zhao
 * Create date: 2011-06-13
 * Contact: zhangmp@bc.edu
 * Last revise: 2011-06-24 
 */

#include "bam.h"
#include "hmm.h"

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
	matrix_array[0][9] = d/L;	/* S -> M_k+1 */
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

float** emission_init (char* ref)
{  
  /* pad A     C     pad G     pad pad pad T     pad pad pad pad pad pad N   
	{0,  1,    0,    0,  0,    0,  0,  0,  0,    0,  0,  0,  0,  0,  0,  0.25},	A 
	{0,  0,    1,    0,  0,    0,  0,  0,  0,    0,  0,  0,  0,  0,  0,  0.25},	C 
	{0,  0,    0,    0,  1,    0,  0,  0,  0,    0,  0,  0,  0,  0,  0,  0.25},	G 
	{0,  0,    0,    0,  0,    0,  0,  0,  1,    0,  0,  0,  0,  0,  0,  0.25},	T 
	{0,  0,    0,    0,  0.5,  0,  0,  0,  0.5,  0,  0,  0,  0,  0,  0,  0.5},	K 
	{0,  0.5,  0.5,  0,  0,    0,  0,  0,  0,    0,  0,  0,  0,  0,  0,  0.5},	M
	{0,  0.5,  0,    0,  0.5,  0,  0,  0,  0,    0,  0,  0,  0,  0,  0,  0.5},	R
	{0,  0,    0.5,  0,  0,    0,  0,  0,  0.5,  0,  0,  0,  0,  0,  0,  0.5},	Y
	{0,  0,    0.5,  0,  0.5,  0,  0,  0,  0,    0,  0,  0,  0,  0,  0,  0.5},	S
	{0,  0,    0.33, 0,  0.33, 0,  0,  0,  0.33, 0,  0,  0,  0,  0,  0,  0.75},	B
	{0,  0.33, 0.33, 0,  0.33, 0,  0,  0,  0,    0,  0,  0,  0,  0,  0,  0.75},	V
	{0,  0.33, 0.33, 0,  0,    0,  0,  0,  0.33, 0,  0,  0,  0,  0,  0,  0.75},	H
	{0,  0.33, 0,    0,  0.33, 0,  0,  0,  0.33, 0,  0,  0,  0,  0,  0,  0.75},	D
	{0,  0.25, 0.25, 0,  0.25, 0,  0,  0,  0.25, 0,  0,  0,  0,  0,  0,  1},	N
	{0,  0,    0,    0,  0,    0,  0,  0,  0,    0,  0,  0,  0,  0,  0,  0},	X
	{0,  0,    0,    0,  0,    0,  0,  0,  0,    0,  0,  0,  0,  0,  0,  0},	pad */
	int32_t i;
	int32_t ref_len = strlen(ref);
	float ** array = (float**)calloc(ref_len, sizeof(float*));
	for (i = 0; i < ref_len; i ++) {
		array[i] = (float*)calloc(16, sizeof(float));
		switch (ref[i]) {
			case 'A':
			case 'a':
				array[i][1] = 1; array[i][15] = 0.25; array[i][2] = array[i][4] = array[i][8] = 0;	
				break;
			case 'C':
			case 'c':
				array[i][2] = 1; array[i][15] = 0.25; array[i][1] = array[i][4] = array[i][8] = 0;	
				break;
			case 'G':
			case 'g':
				array[i][4] = 1; array[i][15] = 0.25; array[i][1] = array[i][2] = array[i][8] = 0;	
				break;
			case 'T':
			case 't':
				array[i][8] = 1; array[i][15] = 0.25; array[i][1] = array[i][2] = array[i][4] = 0;	
				break;
			case 'K': /* G or T */
			case 'k':
				array[i][4] = array[i][8] = array[i][15] = 0.5; array[i][1] = array[i][2] = 0;	
				break;
			case 'M': /* A or C */
			case 'm':
				array[i][1] = array[i][2] = array[i][15] = 0.5; array[i][4] = array[i][8] = 0;	
				break;
			case 'R': /* A or G */
			case 'r':
				array[i][1] = array[i][4] = array[i][15] = 0.5; array[i][2] = array[i][8] = 0;	
				break;
			case 'Y': /* C or T */
			case 'y':
				array[i][2] = array[i][8] = array[i][15] = 0.5; array[i][1] = array[i][4] = 0;	
				break;
			case 'S': /* A or T */
			case 's':
				array[i][2] = array[i][4] = array[i][15] = 0.5; array[i][1] = array[i][8] = 0;	
				break;
			case 'B': /* C or G or T */
			case 'b':
				array[i][2] = array[i][4] = array[i][8] = 0.33; array[i][15] = 0.75; array[i][1] = 0;	
				break;
			case 'V': /* A or C or G */
			case 'v':
				array[i][1] = array[i][2] = array[i][4] = 0.33; array[i][15] = 0.75; array[i][8] = 0;	
				break;
			case 'H': /* A or C or T */
			case 'h':
				array[i][1] = array[i][2] = array[i][8] = 0.33; array[i][15] = 0.75; array[i][4] = 0;	
				break;
			case 'D': /* A or G or T */
			case 'd':
				array[i][1] = array[i][4] = array[i][8] = 0.33; array[i][15] = 0.75; array[i][2] = 0;	
				break;
			case 'N': /* any */
			case 'n':
				array[i][1] = array[i][2] = array[i][4] = array[i][8] = 0.25; array[i][15] = 1;	
				break;
			case 'X': /* any, x mask on Y chromosome */
			case 'x':
				array[i][1] = array[i][2] = array[i][4] = array[i][8] = array[i][15] = 0;	
				break;
			default:
				array[i][1] = array[i][2] = array[i][4] = array[i][8] = array[i][15] = 0;	
				fprintf(stderr, "Wrong reference sequence. \n");
				break;
		}
	}
	return array;
}

void emission_destroy (float** array, const int32_t L)
{
	int32_t i;
	for (i = 0; i < L; i ++) {
		free(array[i]);
	}
	free(array);
}

void forward_backward (float** transition, float** emission, char* ref, char* read)
{
	/* forward algorithm */
	int32_t i;	 /* iter of read */ 
	int32_t k;	 /* iter of reference */
	int32_t ref_len = strlen(ref);
	int32_t read_len = 2*strlen(read); 
	
	/* read: 0-based; reference: 1-based */
	double match[read_len][ref_len + 1];
	double insertion[read_len][ref_len + 1];
	double deletion[read_len][ref_len + 1];
	double s[read_len]; /* scaling factor to avoid underflow */

	double f = 0;

	match[0][0] = deletion[0][0] = deletion[0][1] = 0; /* no M_0, D_0, D_1 states */
	insertion[0][0] = insertion[0][1] = 0.25 * transition[0][10]; /* f_1_I0, f_1_I1 */
	s[0] = insertion[0][0] + insertion[0][1];
	for (k = 1; k <= ref_len; k ++) {
		match[0][k] = emission[k - 1][bam1_seqi(read, 0)] * transition[k - 1][9];
		insertion[0][k] = 0.25 * transition[k][10];
/*		fprintf(stderr, "transition[%d][10]: %f\n", k, transition[k][10]);*/
		s[0] += match[0][k] + insertion[0][k];
	/*	fprintf(stderr, "match[0][%d]: %f\tinsertion[0][%d]: %f\n", k, match[0][k], k, insertion[0][k]); 
*/
	}

	/* rescale */
	for (k = 0; k <= ref_len; k ++) {
		match[0][k] /= s[0];
		insertion[0][k] /= s[0];
	}
	
	for (i = 1; i < read_len; i ++) {
		match[i][0] = deletion[i][0] = deletion[i][1] = 0; /* no M_0, D_0, D_1 states */
		insertion[i][0] = 0.25 * transition[0][5] * insertion[i - 1][0]; /* f_i_I0; i = 2 ~ l */
		
		match[i][1] = emission[0][bam1_seqi(read, i)] * transition[0][4] * 
		insertion[i - 1][0]; /* f_i_M1; i = 1 ~ l */
		
		insertion[i][1] = 0.25 * (transition[1][1] * match[i - 1][1] + transition[1][5] * 
		insertion[i - 1][1]); /* f_i_I1; i = 2 ~ l */
		s[i] = insertion[i][0] + match[i][1] + insertion[i][1];
		for (k = 2; k <= ref_len; k ++) {
			match[i][k] = emission[k - 1][bam1_seqi(read, i)] * (transition[k - 1][0] *
		    match[i - 1][k - 1] + transition[k - 1][4] * insertion[i - 1][k - 1] + transition[k - 1][7] *
			deletion[i - 1][k - 1]);
/*			fprintf(stderr, "**e: %f, t: %f, match[%d][%d]: %f, t: %f, i: %f, t: %f, deletion[%d][%d]: %f**\n", emission[k - 1][bam1_seqi(read, i)], transition[k - 1][0], 
		    i - 1, k - 1, match[i - 1][k - 1], transition[k - 1][4], insertion[i - 1][k - 1], transition[k - 1][7], i - 1, k - 1,  
			deletion[i - 1][k - 1]);
*/
			insertion[i][k] = 0.25 * (transition[k][1] * match[i - 1][k] + transition[k][5] * insertion[i - 1][k]);
/*			fprintf(stderr, "*transition[%d][1]: %f\tmatch[%d][%d]: %f\ttransition[%d][5]: %f\tinsertion[%d][%d]: %f*\n", k, transition[k][1], i - 1, k, match[i - 1][k], k, transition[k][5], i - 1, k, insertion[i - 1][k]);
*/			
			deletion[i][k] = transition[k - 1][2] * match[i][k - 1] + transition[k - 1][8] * deletion[i][k - 1];
/*			fprintf(stderr, "*t: %f, match[%d][%d]: %f, t: %f, insertion[%d][%d]: %f*\n", transition[k - 1][2], i, k - 1, match[i][k - 1], transition[k - 1][8], i, k - 1, insertion[i][k - 1]);
			fprintf(stderr, "match[%d][%d]: %f\tinsertion[%d][%d]: %f\tdeletion[%d][%d]: %f\n\n", i, k, match[i][k], i, k, insertion[i][k], i, k, deletion[i][k]);
*/
			s[i] += match[i][k] + insertion[i][k] + deletion[i][k];
		}
		
		/* rescale */
		for (k = 0; k <= ref_len; k ++) {
			match[i][k] /= s[i];
			insertion[i][k] /= s[i];
			deletion[i][k] /= s[i];
		}

	}

	/* sum of all forward path */
	for (k = 0; k <= ref_len; k ++) {
		f += transition[k][3] * match[read_len - 1][k] + transition[k][6] * insertion[read_len - 1][k];
/*		fprintf(stderr, "match[read_len - 1][%d]: %f\tinsertion[read_len - 1][%d]: %f\n", k, match[read_len - 1][k], k, insertion[read_len - 1][k]); 
*/
	}
	
/*	for (i = 0; i < read_len; i ++) {
		for (k = 0; k <= ref_len; k ++) {
			fprintf(stderr, "insertion[%d][%d]: %f\t", i, k, insertion[i][k]);
		}
		fprintf(stderr, "\n");
	}*/
	fprintf(stderr, "forward: %f\n", f);

	/* backward algorithm */
/*	int32_t i;	  iter of read  
	int32_t k;	  iter of reference
	int32_t ref_len = strlen(ref);
	int32_t read_len = 2*strlen(read);

	 read: 0-based; reference: 1-based 
	double match[read_len][ref_len + 1];
	double insertion[read_len][ref_len + 1];
	double deletion[read_len][ref_len + 1];
*/
	double b = 0;

	match[read_len - 1][0] = 0; /* no M_0 state */
	/*s[read_len - 1] = */insertion[read_len - 1][0] = transition[0][6];
	for (k = 1; k <= ref_len; k ++) {
		match[read_len - 1][k] = transition[k][3];
		insertion[read_len - 1][k] = transition[k][6];
/*		s[read_len - 1] += match[read_len - 1][k] + insertion[read_len - 1][k];*/
	}

	/* rescale */
	for (k = 0; k <= ref_len; k ++) {
		match[read_len - 1][k] /= s[read_len - 1];
		insertion[read_len - 1][k] /= s[read_len - 1];
/*		fprintf(stderr, "match[read_len - 1][k]: %f\tinsertion[read_len - 1][k]: %f\n", match[read_len - 1][k], insertion[read_len - 1][k]);
*/
	}
	
	for (i = read_len - 2; i > 0; i --) {
		match[i][ref_len] = 0.25 * transition[ref_len][1] * insertion[i + 1][ref_len];
		insertion[i][ref_len] = 0.25 * transition[ref_len][5] * insertion[i + 1][ref_len];
		deletion[i][ref_len] = 0;
	/*	s[i] = match[i][ref_len] + insertion[i][ref_len];*/
	 
		for (k = ref_len  - 1; k > 0; k --) {
			match[i][k] = emission[k][bam1_seqi(read, i + 1)] * transition[k][0] * match[i + 1][k + 1] +
			0.25 * transition[k][1] * insertion[i + 1][k] + transition[k][2] * deletion[i][k + 1];

			insertion[i][k] = emission[k][bam1_seqi(read, i + 1)] * transition[k][4] * 
			match[i + 1][k + 1] + 0.25 * transition[k][5] * insertion[i + 1][k];

			deletion[i][k] = emission[k][bam1_seqi(read, i + 1)] * transition[k][7] * 
			match[i + 1][k + 1] + transition[k][8] * deletion[i][k + 1];
/*			fprintf(stderr, "match[%d][%d]: %f, insertion[%d][%d]: %f, deletion[%d][%d]: %f\n", i, k, match[i][k], i, k, insertion[i][k], i, k, deletion[i][k]);
*/
		/*	s[i] += match[i][k] + insertion[i][k] + deletion[i][k];*/
		}
		insertion[i][0] = emission[0][bam1_seqi(read, i + 1)] * transition[0][4] * match[i + 1][1] +
		0.25 * transition[0][5] * insertion[i + 1][0];
		match[i][0] = deletion[i][0] = 0; /* no M_0, D_0 states */
	/*	s[i] += insertion[i][0];*/

		/* rescale */
		for (k = 0; k <= ref_len; k ++) {
			match[i][k] /= s[i];
			insertion[i][k] /= s[i];
			deletion[i][k] /= s[i];
		}

	}

/*	s[0] = 0;*/
	deletion[0][ref_len] = 0;
	match[0][ref_len] = 0;
	insertion[0][ref_len] = 0.25 * transition[ref_len][5] * insertion[1][ref_len];
	for (k = ref_len - 1; k >= 0; k --) {
		match[0][k] = emission[k][bam1_seqi(read, 1)] * transition[k][0] * match[1][k + 1] +
		0.25 * transition[k][1] * insertion[i + 1][k] + transition[k][2] * deletion[0][k + 1];

		insertion[0][k] = emission[k][bam1_seqi(read, 1)] * transition[k][4] * 
		match[1][k + 1] + 0.25 * transition[k][5] * insertion[1][k];

		deletion[0][k] = 0;
/*		fprintf(stderr, "match[0][%d]: %f\tinsertion[0][%d]: %f\n", k, match[0][k], k, insertion[0][k]);*/
	/*	s[0] += match[0][k] + insertion[0][k];*/
		/* sum of all backward path */
	}

	/* rescale */
	for (k = 1; k <= ref_len; k ++) {
		match[0][k] /= s[0];
		insertion[0][k] /= s[0];

		b += emission[k - 1][bam1_seqi(read, 0)] * transition[k - 1][9] * match[0][k] + 0.25 *
		transition[k][10] * insertion[0][k];
/*		fprintf(stderr, "match[0][%d]: %f\tinsertion[0][%d]: %f\n", k, match[0][k], k, insertion[0][k]);*/
	}
	insertion[0][0] /= s[0];

	b += 0.25 * transition[0][10] * insertion[0][0]; 	
/*	fprintf(stderr, "insertion[0][0]: %f\n", insertion[0][0]);
*/
	fprintf(stderr, "backward: %f\n", b);
} 
