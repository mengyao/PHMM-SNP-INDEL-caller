/*
 * hmm.c: Banded semi-global HMM
 * Author: Mengyao Zhao
 * Create date: 2011-06-13
 * Contact: zhangmp@bc.edu
 * Last revise: 2011-07-11 
 */

#include <math.h>
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
	double f_match[read_len][ref_len + 1];
	double f_insertion[read_len][ref_len + 1];
	double f_deletion[read_len][ref_len + 1];
	double s[read_len + 1]; /* scaling factor to avoid underflow */

	double f = 0;

	f_match[0][0] = f_deletion[0][0] = f_deletion[0][1] = 0; /* no M_0, D_0, D_1 states */
	f_insertion[0][0] = f_insertion[0][1] = 0.25 * transition[0][10]; /* f_1_I0, f_1_I1 */
	s[0] = f_insertion[0][0] + f_insertion[0][1];
	for (k = 1; k <= ref_len; k ++) {
		f_match[0][k] = emission[k - 1][bam1_seqi(read, 0)] * transition[k - 1][9];
		f_insertion[0][k] = 0.25 * transition[k][10];
		s[0] += f_match[0][k] + f_insertion[0][k];
	}

	/* rescale */
	for (k = 0; k <= ref_len; k ++) {
		f_match[0][k] /= s[0];
		f_insertion[0][k] /= s[0];
	}
	
	for (i = 1; i < read_len; i ++) {
		f_match[i][0] = f_deletion[i][0] = f_deletion[i][1] = 0; /* no M_0, D_0, D_1 states */
		f_insertion[i][0] = 0.25 * transition[0][5] * f_insertion[i - 1][0]; /* f_i_I0; i = 2 ~ l */
		
		f_match[i][1] = emission[0][bam1_seqi(read, i)] * transition[0][4] * 
		f_insertion[i - 1][0]; /* f_i_M1; i = 1 ~ l */
		
		f_insertion[i][1] = 0.25 * (transition[1][1] * f_match[i - 1][1] + transition[1][5] * 
		f_insertion[i - 1][1]); /* f_i_I1; i = 2 ~ l */
		s[i] = f_insertion[i][0] + f_match[i][1] + f_insertion[i][1];
		for (k = 2; k <= ref_len; k ++) {
			f_match[i][k] = emission[k - 1][bam1_seqi(read, i)] * (transition[k - 1][0] *
		    f_match[i - 1][k - 1] + transition[k - 1][4] * f_insertion[i - 1][k - 1] + transition[k - 1][7] *
			f_deletion[i - 1][k - 1]);
			
			f_insertion[i][k] = 0.25 * (transition[k][1] * f_match[i - 1][k] + transition[k][5] * f_insertion[i - 1][k]);
			
			f_deletion[i][k] = transition[k - 1][2] * f_match[i][k - 1] + transition[k - 1][8] * f_deletion[i][k - 1];
			
			s[i] += f_match[i][k] + f_insertion[i][k] + f_deletion[i][k];
		}
		
		/* rescale */
		for (k = 0; k <= ref_len; k ++) {
			f_match[i][k] /= s[i];
			f_insertion[i][k] /= s[i];
			f_deletion[i][k] /= s[i];
		}
	}

	/* sum of all forward path */
	for (k = 0; k <= ref_len; k ++) {
		f += transition[k][3] * f_match[read_len - 1][k] + transition[k][6] * f_insertion[read_len - 1][k];
	}
	
	s[read_len] = f;
	f /= s[read_len];
	fprintf(stderr, "forward: %f\n", f);

	{// compute the log likelihood
		double p = 1., Pr1 = 0., Pr;
		for (i = 0; i <= read_len; ++i) {
			p *= s[i];
			if (p < 1e-100) Pr += -4.343 * log(p), p = 1.;
		}
		Pr1 += -4.343 * log(p * (ref_len + 1) * read_len);
		Pr = (int)(Pr1 + .499);
		fprintf (stderr, "log likelihood: %f\n", Pr);
	}

	/* backward algorithm */
	/* read: 0-based; reference: 1-based */
	double b_match[read_len][ref_len + 1];
	double b_insertion[read_len][ref_len + 1];
	double b_deletion[read_len][ref_len + 1];
	double b = 0;

	b_match[read_len - 1][0] = 0; /* no M_0 state */
	b_insertion[read_len - 1][0] = transition[0][6];
	for (k = 1; k <= ref_len; k ++) {
		b_match[read_len - 1][k] = transition[k][3];
		b_insertion[read_len - 1][k] = transition[k][6];
	}

	/* rescale */
	for (k = 0; k <= ref_len; k ++) {
		/* b_0_E needs to be rescaled by s[read_len] */
		b_match[read_len - 1][k] = b_match[read_len - 1][k] / s[read_len - 1] / s[read_len];
		b_insertion[read_len - 1][k] = b_insertion[read_len - 1][k] / s[read_len - 1] / s[read_len];
	}
	
	for (i = read_len - 2; i > 0; i --) {
		b_match[i][ref_len] = 0.25 * transition[ref_len][1] * b_insertion[i + 1][ref_len];

		b_insertion[i][ref_len] = 0.25 * transition[ref_len][5] * b_insertion[i + 1][ref_len];

		b_deletion[i][ref_len] = 0;
	 
		for (k = ref_len  - 1; k > 0; k --) {
			b_match[i][k] = emission[k][bam1_seqi(read, i + 1)] * transition[k][0] * b_match[i + 1][k + 1] +
			0.25 * transition[k][1] * b_insertion[i + 1][k] + transition[k][2] * b_deletion[i][k + 1];

			b_insertion[i][k] = emission[k][bam1_seqi(read, i + 1)] * transition[k][4] * 
			b_match[i + 1][k + 1] + 0.25 * transition[k][5] * b_insertion[i + 1][k];

			b_deletion[i][k] = emission[k][bam1_seqi(read, i + 1)] * transition[k][7] * 
			b_match[i + 1][k + 1] + transition[k][8] * b_deletion[i][k + 1];
		}
		b_insertion[i][0] = emission[0][bam1_seqi(read, i + 1)] * transition[0][4] * b_match[i + 1][1] +
		0.25 * transition[0][5] * b_insertion[i + 1][0];
		b_match[i][0] = b_deletion[i][0] = 0; /* no M_0, D_0 states */

		/* rescale */
		for (k = 0; k <= ref_len; k ++) {
			b_match[i][k] /= s[i];
			b_insertion[i][k] /= s[i];
			b_deletion[i][k] /= s[i];
		}
	}

	b_deletion[0][ref_len] = 0;
	b_match[0][ref_len] = 0;
	b_insertion[0][ref_len] = 0.25 * transition[ref_len][5] * b_insertion[1][ref_len];
	for (k = ref_len - 1; k >= 0; k --) {
		b_match[0][k] = emission[k][bam1_seqi(read, 1)] * transition[k][0] * b_match[1][k + 1] +
		0.25 * transition[k][1] * b_insertion[i + 1][k] + transition[k][2] * b_deletion[0][k + 1];

		b_insertion[0][k] = emission[k][bam1_seqi(read, 1)] * transition[k][4] * 
		b_match[1][k + 1] + 0.25 * transition[k][5] * b_insertion[1][k];

		b_deletion[0][k] = 0;
	}

	/* rescale */
	for (k = 1; k <= ref_len; k ++) {
		b_match[0][k] /= s[0];
		b_insertion[0][k] /= s[0];

		b += emission[k - 1][bam1_seqi(read, 0)] * transition[k - 1][9] * b_match[0][k] + 0.25 *
		transition[k][10] * b_insertion[0][k];
	}
	b_insertion[0][0] /= s[0];

	b += 0.25 * transition[0][10] * b_insertion[0][0]; 	
	fprintf(stderr, "backward: %f\n", b);
}

void baum_welch (char* ref_seq, int len, bamFile fp, char* bai, int32_t target_num, int32_t begin, int32_t end) /* 0-based coordinate */ 
{
	bam1_t* b = bam_init1();
	bam_index_t* idx = bam_index_load(bai);
	bam_iter_t bam_iter = bam_iter_query(idx, target_num, begin, end);
	fprintf (stdout, "reference sequence: %s\n", ref_seq); 
	float** matrix_array = transition_init (0.3, 0.5, 0.2, 0.5, 0.5, len);
	float** emission = emission_init(ref_seq);
	
	while (bam_iter_read (fp, bam_iter, b) >= 0) {
		char* read_seq = (char*)bam1_seq(b);
		char* read_name = bam1_qname(b);
		fprintf (stdout, "read name: %s\n", read_name);
		forward_backward (matrix_array, emission, ref_seq, read_seq);
	}
	bam_iter_destroy(bam_iter);
	emission_destroy(emission, len);
	transition_destroy(matrix_array, len);
	bam_index_destroy(idx);
	bam_destroy1(b);
} 
