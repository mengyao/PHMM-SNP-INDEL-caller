/*
 * hmm.c: Semi-global HMM
 * Author: Mengyao Zhao
 * Create date: 2011-06-13
 * Contact: zhangmp@bc.edu
 * Last revise: 2011-08-07 
 */

#include <math.h>
#include "bam.h"
#include "hmm.h"

double** transition_init (const double a, const double b, const double r, const double c, const double d, const int32_t L)
{
	double** matrix_array = calloc (L + 1, sizeof(double*));
	int32_t i;
	for (i = 0; i < L + 1; i ++) {
		matrix_array[i] = calloc (16, sizeof(double));
	}

	/* k = 0: inseart before the reference */	
/*	for (i = 0; i < 16; i ++) {
		matrix_array[0][i] = 0;
	}*/
	matrix_array[0][4] = (1 - c)*r;	/* I_k -> M_k+1 */
	matrix_array[0][5] = c*r;	/* I_k -> I_k */
	matrix_array[0][6] = r;	/* I_k -> E */
	matrix_array[0][9] = d/L;	/* S -> M_k+1 */
	matrix_array[0][10] = (1 - d)/(L + 1);	/* S -> I_k */

	/* k = 1 ... L - 2 */
	for (i = 1; i < L - 1; i ++) {
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
	
		/* Sum of the out edges of S equals to 1 */	
		matrix_array[i][9] = d/L;	/* S -> M_k+1 */
		matrix_array[i][10] = (1 - d)/(L + 1);	/* S -> I_k */
	}	
	/* State D_1 doesn't exist. */
	matrix_array[1][7] = 0;
	matrix_array[1][8] = 0;	

	/* k = L - 1 */
	/* Sum of the following 4 lines equals to 1. */
	matrix_array[L - 1][0] = (1 - a)*(1 - r);	/* M_k -> M_k+1 */
	matrix_array[L - 1][1] = a*(1 - r);	/* M_k -> I_k */
	matrix_array[L - 1][2] = 0;	/* M_k -> D_k+1, no D_L state */
	matrix_array[L - 1][3] = r;	/* M_k -> E */

	/* Sum of the following 3 lines equals to 1. */
	matrix_array[L - 1][4] = (1 - c)*r;	/* I_k -> M_k+1 */
	matrix_array[L - 1][5] = c*r;	/* I_k -> I_k */
	matrix_array[L - 1][6] = r;	/* I_k -> E */

	/* Sum of the following 2 lines equals to 1. */
	matrix_array[L - 1][7] = 1;	/* D_k -> M_k+1 */
	matrix_array[L - 1][8] = 0;		/* D_k -> D_k+1, no D_L state */

	/* Sum of the out edges of S equals to 1 */	
	matrix_array[L - 1][9] = d/L;	/* S -> M_k+1 */
	matrix_array[L - 1][10] = (1 - d)/(L + 1);	/* S -> I_k */
	
	/* k = L */
/*	for (i = 0; i < 16; i ++) {
		matrix_array[L][i] = 0;
	}*/
	matrix_array[L][1] = a*(1 - r);	/* M_k -> I_k */
	matrix_array[L][3] = r;	/* M_k -> E */
	matrix_array[L][5] = 1 - r;	/* I_k -> I_k */
	matrix_array[L][6] = r;	/* I_k -> E */
	matrix_array[L][10] = (1 - d)/(L + 1);	/* S -> I_k */
	
	return matrix_array;
} 

void transition_destroy (double** matrix_array, const int32_t L)
{
	int32_t i;
	for (i = 0; i < L + 1; i ++) {
		free(matrix_array[i]);
	}
	free(matrix_array);
}

double** emission_init (char* ref)
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
	double ** array = (double**)calloc(ref_len, sizeof(double*));
	for (i = 0; i < ref_len; i ++) {
		array[i] = (double*)calloc(16, sizeof(double));
		switch (ref[i]) {
			case 'A':
			case 'a':
				array[i][1] = 1; array[i][15] = 0.25; array[i][2] = array[i][4] = array[i][8] = 0.001;	
				break;
			case 'C':
			case 'c':
				array[i][2] = 1; array[i][15] = 0.25; array[i][1] = array[i][4] = array[i][8] = 0.001;	
				break;
			case 'G':
			case 'g':
				array[i][4] = 1; array[i][15] = 0.25; array[i][1] = array[i][2] = array[i][8] = 0.001;	
				break;
			case 'T':
			case 't':
				array[i][8] = 1; array[i][15] = 0.25; array[i][1] = array[i][2] = array[i][4] = 0.001;	
				break;
			case 'N': /* any */
			case 'n':
				array[i][1] = array[i][2] = array[i][4] = array[i][8] = 0.25; array[i][15] = 1;	
				break;
			case 'X': /* any, x mask on Y chromosome */
			case 'x':
				array[i][1] = array[i][2] = array[i][4] = array[i][8] = array[i][15] = 0.001;	
				break;
			case 'K': /* G or T */
			case 'k':
				array[i][4] = array[i][8] = array[i][15] = 0.5; array[i][1] = array[i][2] = 0.001;	
				break;
			case 'M': /* A or C */
			case 'm':
				array[i][1] = array[i][2] = array[i][15] = 0.5; array[i][4] = array[i][8] = 0.001;	
				break;
			case 'R': /* A or G */
			case 'r':
				array[i][1] = array[i][4] = array[i][15] = 0.5; array[i][2] = array[i][8] = 0.001;	
				break;
			case 'Y': /* C or T */
			case 'y':
				array[i][2] = array[i][8] = array[i][15] = 0.5; array[i][1] = array[i][4] = 0.001;	
				break;
			case 'S': /* A or T */
			case 's':
				array[i][2] = array[i][4] = array[i][15] = 0.5; array[i][1] = array[i][8] = 0.001;	
				break;
			case 'B': /* C or G or T */
			case 'b':
				array[i][2] = array[i][4] = array[i][8] = 0.33; array[i][15] = 0.75; array[i][1] = 0.001;	
				break;
			case 'V': /* A or C or G */
			case 'v':
				array[i][1] = array[i][2] = array[i][4] = 0.33; array[i][15] = 0.75; array[i][8] = 0.001;	
				break;
			case 'H': /* A or C or T */
			case 'h':
				array[i][1] = array[i][2] = array[i][8] = 0.33; array[i][15] = 0.75; array[i][4] = 0.001;	
				break;
			case 'D': /* A or G or T */
			case 'd':
				array[i][1] = array[i][4] = array[i][8] = 0.33; array[i][15] = 0.75; array[i][2] = 0.001;	
				break;
			default:
				/*array[i][1] = array[i][2] = array[i][4] = array[i][8] = array[i][15] = 0;*/	
				fprintf(stderr, "Wrong reference sequence. \n");
				exit (1);
				break;
		}
	}
	return array;
}

void emission_destroy (double** array, const int32_t L)
{
	int32_t i;
	for (i = 0; i < L; i ++) {
		free(array[i]);
	}
	free(array);
}

double forward_backward (double** transition, double** emission, char* ref, uint8_t* read, int32_t read_len, fb* f, fb* b, double* s)
{
	/*-------------------*
	 * forward algorithm *
	 *-------------------*/
	int32_t i;	 /* iter of read */ 
	int32_t k;	 /* iter of reference */
	int32_t ref_len = strlen(ref);
	/* double pp = 0;	 Debug: posterior probability */

	/* f_0_S = s_begin = 1 */	
	f->match[0][0] = f->deletion[0][0] = f->deletion[0][1] = 0; /* no M_0, D_0, D_1 states */
	f->insertion[0][0] = f->insertion[0][1] = 0.25 * transition[0][10]; /* f_1_I0, f_1_I1 */
	s[0] = f->insertion[0][0] + f->insertion[0][1];
	for (k = 1; k <= ref_len; k ++) {
		f->match[0][k] = emission[k - 1][bam1_seqi(read, 0)] * transition[k - 1][9];
		f->insertion[0][k] = 0.25 * transition[k][10];
		s[0] += f->match[0][k] + f->insertion[0][k];
	}

	/* rescale */
	for (k = 0; k <= ref_len; k ++) {
		f->match[0][k] /= s[0];
		f->insertion[0][k] /= s[0];
	}
	
	for (i = 1; i < read_len; i ++) {
		f->match[i][0] = f->deletion[i][0] = f->deletion[i][1] = 0; /* no M_0, D_0, D_1 states */
		f->insertion[i][0] = 0.25 * transition[0][5] * f->insertion[i - 1][0]; /* f_i_I0; i = 2 ~ l */
		
		f->match[i][1] = emission[0][bam1_seqi(read, i)] * transition[0][4] * 
		f->insertion[i - 1][0]; /* f_i_M1; i = 1 ~ l */
		
		f->insertion[i][1] = 0.25 * (transition[1][1] * f->match[i - 1][1] + transition[1][5] * 
		f->insertion[i - 1][1]); /* f_i_I1; i = 2 ~ l */
		s[i] = f->insertion[i][0] + f->match[i][1] + f->insertion[i][1];
		
		for (k = 2; k < ref_len; k ++) {
			f->match[i][k] = emission[k - 1][bam1_seqi(read, i)] * (transition[k - 1][0] *
		    f->match[i - 1][k - 1] + transition[k - 1][4] * f->insertion[i - 1][k - 1] + transition[k - 1][7] *
			f->deletion[i - 1][k - 1]);
			
			f->insertion[i][k] = 0.25 * (transition[k][1] * f->match[i - 1][k] + transition[k][5] * f->insertion[i - 1][k]);
			
			f->deletion[i][k] = transition[k - 1][2] * f->match[i][k - 1] + transition[k - 1][8] * f->deletion[i][k - 1];
			
			s[i] += f->match[i][k] + f->insertion[i][k] + f->deletion[i][k];
		}

		f->match[i][ref_len] = emission[ref_len - 1][bam1_seqi(read, i)] * (transition[ref_len - 1][0] *
		f->match[i - 1][ref_len - 1] + transition[ref_len - 1][4] * f->insertion[i - 1][ref_len - 1] + transition[ref_len - 1][7] *
		f->deletion[i - 1][ref_len - 1]);
		
		f->insertion[i][ref_len] = 0.25 * (transition[ref_len][1] * f->match[i - 1][ref_len] + transition[ref_len][5] * f->insertion[i - 1][ref_len]);
		
		f->deletion[i][ref_len] = 0;	/* no D_L state */
		
		s[i] += f->match[i][ref_len] + f->insertion[i][ref_len];
		
		/* rescale */
		for (k = 0; k <= ref_len; k ++) {
			f->match[i][k] /= s[i];
			f->insertion[i][k] /= s[i];
			f->deletion[i][k] /= s[i];
		}
	}

	/* sum of all forward path */
	f->final = 0;
	for (k = 0; k <= ref_len; k ++) {
		f->final += transition[k][3] * f->match[read_len - 1][k] + transition[k][6] * f->insertion[read_len - 1][k];
	}
	
	s[read_len] = f->final;
	f->final /= s[read_len];

	/*--------------------*
	 * backword algorithm *
	 *--------------------*/
	b->match[read_len - 1][0] = 0; /* no M_0 state */
	b->insertion[read_len - 1][0] = transition[0][6] / s[read_len];
	for (k = 1; k <= ref_len; k ++) {
		b->match[read_len - 1][k] = transition[k][3] / s[read_len];
		b->insertion[read_len - 1][k] = transition[k][6] / s[read_len];
	}

	/* rescale */
	if (read_len > 1) {	
		for (k = 0; k <= ref_len; k ++) {
			/* b_0_E needs to be rescaled by s[read_len] */
			b->match[read_len - 1][k] = b->match[read_len - 1][k] / s[read_len - 1] /* s[read_len] */;
			b->insertion[read_len - 1][k] = b->insertion[read_len - 1][k] / s[read_len - 1] /* s[read_len] */;
			/* Debug: posterior probability */
			/* pp += b->match[read_len - 1][k] * f->match[read_len - 1][k] + b->insertion[read_len - 1][k] * f->insertion[read_len - 1][k]; */
		}
	/*	pp *= s[read_len - 1];
		fprintf (stderr, "pp: %f\n", pp);*/
	
		for (i = read_len - 2; i > 0; i --) {
			b->match[i][ref_len] = 0.25 * transition[ref_len][1] * b->insertion[i + 1][ref_len];

			b->insertion[i][ref_len] = 0.25 * transition[ref_len][5] * b->insertion[i + 1][ref_len];

			b->deletion[i][ref_len] = 0;
		 
			b->match[i][ref_len - 1] = emission[ref_len - 1][bam1_seqi(read, i + 1)] * transition[ref_len - 1][0] * b->match[i + 1][ref_len] +
			0.25 * transition[ref_len - 1][1] * b->insertion[i + 1][ref_len - 1];	/* no D_L state */

			b->insertion[i][ref_len - 1] = emission[ref_len - 1][bam1_seqi(read, i + 1)] * transition[ref_len - 1][4] * 
			b->match[i + 1][ref_len] + 0.25 * transition[ref_len - 1][5] * b->insertion[i + 1][ref_len - 1];

			b->deletion[i][ref_len - 1] = emission[ref_len - 1][bam1_seqi(read, i + 1)] * transition[ref_len - 1][7] * 
			b->match[i + 1][ref_len];	/* no D_L state */
			
			for (k = ref_len  - 2; k > 0; k --) {
				b->match[i][k] = emission[k][bam1_seqi(read, i + 1)] * transition[k][0] * b->match[i + 1][k + 1] +
				0.25 * transition[k][1] * b->insertion[i + 1][k] + transition[k][2] * b->deletion[i][k + 1];

				b->insertion[i][k] = emission[k][bam1_seqi(read, i + 1)] * transition[k][4] * 
				b->match[i + 1][k + 1] + 0.25 * transition[k][5] * b->insertion[i + 1][k];

				b->deletion[i][k] = emission[k][bam1_seqi(read, i + 1)] * transition[k][7] * 
				b->match[i + 1][k + 1] + transition[k][8] * b->deletion[i][k + 1];
			}
			b->deletion[i][1] = 0;	/* no D_1 state */
			b->insertion[i][0] = emission[0][bam1_seqi(read, i + 1)] * transition[0][4] * b->match[i + 1][1] +
			0.25 * transition[0][5] * b->insertion[i + 1][0];
			b->match[i][0] = b->deletion[i][0] = 0; /* no M_0, D_0 states */

			/* rescale */
	/* Debug:	pp = 0; */
			for (k = 0; k <= ref_len; k ++) {
				b->match[i][k] /= s[i];
				b->insertion[i][k] /= s[i];
				b->deletion[i][k] /= s[i];
	/* Debug:	pp += b->match[i][k] * f->match[i][k] + b->insertion[i][k] * f->insertion[i][k]; */
			}
	/* Debug:	pp *= s[i];
			fprintf (stderr, "pp: %f\n", pp);*/
		}

		b->deletion[0][ref_len] = 0;
		b->match[0][ref_len] = 0.25 * transition[ref_len][1] * b->insertion[1][ref_len];
		b->insertion[0][ref_len] = 0.25 * transition[ref_len][5] * b->insertion[1][ref_len];

		b->match[0][ref_len - 1] = emission[ref_len - 1][bam1_seqi(read, 1)] * transition[ref_len - 1][0] * b->match[1][ref_len] +
		0.25 * transition[ref_len - 1][1] * b->insertion[1][ref_len - 1];	/* no D_L state */

		b->insertion[0][ref_len - 1] = emission[ref_len - 1][bam1_seqi(read, 1)] * transition[ref_len - 1][4] * 
		b->match[1][ref_len] + 0.25 * transition[ref_len - 1][5] * b->insertion[1][ref_len - 1];

		b->deletion[0][ref_len - 1] = 0;

		for (k = ref_len - 2; k >= 0; k --) {
			b->match[0][k] = emission[k][bam1_seqi(read, 1)] * transition[k][0] * b->match[1][k + 1] +
			0.25 * transition[k][1] * b->insertion[1][k] + transition[k][2] * b->deletion[0][k + 1];

			b->insertion[0][k] = emission[k][bam1_seqi(read, 1)] * transition[k][4] * 
			b->match[1][k + 1] + 0.25 * transition[k][5] * b->insertion[1][k];

			b->deletion[0][k] = 0;
		}

	}
	/* rescale */
	/* Debug:	pp = 0; */
	b->final = 0;
	for (k = 1; k <= ref_len; k ++) {
		b->match[0][k] /= s[0];
		b->insertion[0][k] /= s[0];
		/* Debug:	pp += b->match[0][k] * f->match[0][k] + b->insertion[0][k] * f->insertion[0][k]; */

		b->final += emission[k - 1][bam1_seqi(read, 0)] * transition[k - 1][9] * b->match[0][k] + 0.25 *
		transition[k][10] * b->insertion[0][k];
	}
	b->insertion[0][0] /= s[0];
	/* Debug:	pp += b->insertion[0][0] * f->insertion[0][0];
	pp *= s[0]; */
	b->final += 0.25 * transition[0][10] * b->insertion[0][0];
	/* Debug:	fprintf (stderr, "pp: %f\nb->final: %g\n", pp, b->final); */ 

	/* Debug: posterior probability for transition 
	{
		double pp_t = 0;
		for (i = 0; i < read_len - 1; i ++) {
			pp_t = 0;
			for (k = 2; k < ref_len; k ++) {
				pp_t += f->deletion[i][k] * transition[k][7] * emission[k][bam1_seqi(read, i + 1)] * b->match[i + 1][k + 1];
			}
			for (k = 1; k < ref_len; k ++) {
				pp_t += f->match[i][k] * transition[k][0] * emission[k][bam1_seqi(read, i + 1)] * b->match[i + 1][k + 1];
				pp_t += f->insertion[i][k] * transition[k][4] * emission[k][bam1_seqi(read, i + 1)] * b->match[i + 1][k + 1];
				pp_t += 0.25 * f->match[i][k] * transition[k][1] * b->insertion[i + 1][k];
				pp_t += 0.25 * f->insertion[i][k] * transition[k][5] * b->insertion[i + 1][k];
			}
			pp_t += f->insertion[i][0] * transition[0][4] * emission[0][bam1_seqi(read, i + 1)] * b->match[i + 1][1];
			pp_t += 0.25 * f->match[i][ref_len] * transition[ref_len][1] * b->insertion[i + 1][ref_len]; 
			pp_t += 0.25 * f->insertion[i][0] * transition[0][5] * b->insertion[i + 1][0];
			pp_t += 0.25 * f->insertion[i][ref_len] * transition[ref_len][5] * b->insertion[i + 1][ref_len];
			fprintf (stderr, "pp_t: %g\n", pp_t);
		}
	} */
		
	{// compute the log likelihood
		double p = 1, Pr1 = 0;
		for (i = 0; i <= read_len; ++i) {
			p *= s[i];
			if (p < 1e-100) Pr1 += log(p), p = 1.;
		}
		Pr1 = log (p);
		return Pr1;
	}
}

void baum_welch (char* ref_seq, int32_t ref_len, reads* r, double df) /* 0-based coordinate */ 
{
	fprintf (stdout, "reference sequence: %s\n", ref_seq); 
	double** transition = transition_init (0.002, 0.98, 0.00067, 0.02, 0.998, ref_len);
	double** emission = emission_init(ref_seq);
	double Pr = 10e100, diff = 1;
	int32_t i, k, j, count = 0;
	double** t = calloc (ref_len + 1, sizeof(double*));
	for (i = 0; i <= ref_len; i ++) {
		t[i] = calloc (16, sizeof(double));
	}		
	double** e = calloc (ref_len, sizeof(double*));
	for (i = 0; i < ref_len; i ++) {
		e[i] = calloc (16, sizeof(double));
	}
	while (diff > df && count < 10) {
		double es = 0; 
		double p = 0;	/* likelihood */
		if (count > 0) {
			for (k = 0; k <= ref_len; k ++) {
				transition[k][0] = t[k][0];
				transition[k][1] = t[k][1]; 
				transition[k][2] = t[k][2];
				transition[k][4] = t[k][4];
				transition[k][5] = t[k][5];
				transition[k][7] = t[k][7];
				transition[k][8] = t[k][8];
			}
			for (k = 0; k < ref_len; k ++) {
				for (i = 0; i < 16; i ++) {
					emission[k][i] = e[k][i];
				}
			}
		}
		/* Initialize new transition and emission matrixes. */
		for (k = 0; k <= ref_len; k ++) {
			t[k][3] = transition[k][3];
			t[k][6] = transition[k][6];
			t[k][9] = transition[k][9];
			t[k][10] = transition[k][10];
			t[k][0] = t[k][1] = t[k][2] = t[k][4] = t[k][5] = t[k][7] = t[k][8] = 0;
		}
		for (k = 0; k < ref_len; k ++) {
			for (i = 0; i < 16; i ++) {
				e[k][i] = 0;
			}
		}
		double s_t[ref_len + 1][3], s_e[ref_len]/*, s_tS = 0*/;
		int32_t total_hl = 0;

		/* Transition and emission matrixes training by a block of reads. */
		for (j = 0; j < r->count; j ++) {
			uint8_t* read_seq = &r->seqs[total_hl];
			total_hl += r->seq_l[j]/2 + r->seq_l[j]%2;
			int32_t read_len = r->seq_l[j];
/*			fprintf (stderr, "read length: %d\nread_seq: ", read_len);
			for (i = 0; i < read_len; i ++) {
				fprintf (stderr, "%d ", bam1_seqi(read_seq, i));
			}
			fprintf (stderr, "\n");
*/
			fb* f = (fb*)calloc(1, sizeof(fb));
			fb* b = (fb*)calloc(1, sizeof(fb));
			f->match = (double**)calloc(read_len, sizeof(double*));
			f->insertion = (double**)calloc(read_len, sizeof(double*));
			f->deletion = (double**)calloc(read_len, sizeof(double*));
			b->match = (double**)calloc(read_len, sizeof(double*));
			b->insertion = (double**)calloc(read_len, sizeof(double*));
			b->deletion = (double**)calloc(read_len, sizeof(double*));
			for (i = 0; i < read_len; i ++) {
				f->match[i] = (double*)calloc(ref_len + 1, sizeof(double));
				f->insertion[i] = (double*)calloc(ref_len + 1, sizeof(double));
				f->deletion[i] = (double*)calloc(ref_len + 1, sizeof(double));
				b->match[i] = (double*)calloc(ref_len + 1, sizeof(double));
				b->insertion[i] = (double*)calloc(ref_len + 1, sizeof(double));
				b->deletion[i] = (double*)calloc(ref_len + 1, sizeof(double));
			}	
			double* s = (double*)calloc(read_len + 1, sizeof(double));			

			p += forward_backward (transition, emission, ref_seq, read_seq, read_len, f, b, s);
			for (k = 0; k < ref_len; k ++) {
				for (i = 0; i < read_len - 1; i ++) {
					t[k][0] += f->match[i][k] * transition[k][0] * emission[k][bam1_seqi(read_seq, i + 1)] 
					* b->match[i + 1][k + 1];	/* M_k -> M_k+1 */

					t[k][1] += 0.25 * f->match[i][k] * transition[k][1] * b->insertion[i + 1][k]; /* M_k -> I_k */

					t[k][2] += f->match[i][k] * transition[k][2] * b->deletion[i][k + 1] * s[i];	/* M_k -> D_k+1 */
					
					t[k][4] += f->insertion[i][k] * transition[k][4] * emission[k][bam1_seqi(read_seq, i + 1)] 
					* b->match[i + 1][k + 1];	/* I_k -> M_k+1 */
			
					t[k][5] += 0.25 * f->insertion[i][k] * transition[k][5] * b->insertion[i + 1][k];	/* I_k -> I_k */

					t[k][7] += f->deletion[i][k] * transition[k][7] * emission[k][bam1_seqi(read_seq, i + 1)] 
					* b->match[i + 1][k + 1];	/* D_k -> M_k+1 */
		
					t[k][8] += f->deletion[i][k] * transition[k][8] * b->deletion[i][k + 1] * s[i];	/* D_k -> D_k+1 */
				}
				/* i = read_len - 1 */
				t[k][2] += f->match[read_len - 1][k] * transition[k][2] * b->deletion[read_len - 1][k + 1] * s[read_len - 1];	/* M_k -> D_k+1 */
				
				t[k][8] += f->deletion[read_len - 1][k] * transition[k][8] * b->deletion[read_len - 1][k + 1] * s[read_len - 1];	/* D_k -> D_k+1 */
				
				for (i = 0; i < read_len; i ++) {
					e[k][bam1_seqi(read_seq, i)] += f->match[i][k + 1] * b->match[i][k + 1] * s[i];
				}
			}
			for (i = 0; i < read_len - 1; i ++) {
				 /* M_k -> I_k */
				t[ref_len][1] += 0.25 * f->match[i][ref_len] * transition[ref_len][1] * b->insertion[i + 1][ref_len];
					
				/* I_k -> I_k */
				t[ref_len][5] += 0.25 * f->insertion[i][ref_len] * transition[ref_len][5] * b->insertion[i + 1][ref_len];
			}
			free(s);
			for (i = 0; i < read_len; i ++) {
				free(f->match[i]);
				free(f->insertion[i]);
				free(f->deletion[i]);
				free(b->match[i]);
				free(b->insertion[i]);
				free(b->deletion[i]);
			}	
			free(f);
			free(b);
		}
		/* Loop ending: Transition and emission matrixes training by a block of reads. */

		/* Set the t with doesn't exist edges to 0 */
		t[0][0] = t[0][1] = t[0][2] = t[0][3] = t[0][7] = t[0][8] = 
		t[1][7] = t[1][8] = t[ref_len - 1][2] = t[ref_len - 1][8] = 
		t[ref_len][0] = t[ref_len][2] = t[ref_len][4] = 
		t[ref_len][7] = t[ref_len][8] = t[ref_len][9] = 0;	

		/* Estimate transition probabilities. */
		s_t[0][1] = t[0][4] + t[0][5] + t[0][6];
		if (s_t[0][1] > 0) {
			t[0][4] /= s_t[0][1];
			t[0][5] /= s_t[0][1];
			t[0][6] /= s_t[0][1];
		}

		for (k = 1; k < ref_len; k ++) {
			s_t[k][0] = t[k][0] + t[k][1] + t[k][2] + t[k][3];
			if (s_t[k][0] > 0) {
				t[k][0] /= s_t[k][0];
				t[k][1] /= s_t[k][0];
				t[k][2] /= s_t[k][0];
				t[k][3] /= s_t[k][0];
			}

			s_t[k][1] = t[k][4] + t[k][5] + t[k][6];
			if (s_t[k][1] > 0) {
				t[k][4] /= s_t[k][1];
				t[k][5] /= s_t[k][1];
				t[k][6] /= s_t[k][1];
			}
	
			s_t[k][2] = t[k][7] + t[k][8];
			if (s_t[k][2] > 0) {
				t[k][7] /= s_t[k][2];
				t[k][8] /= s_t[k][2];
			}
		}

		s_t[ref_len][0] = t[ref_len][1] + t[ref_len][3];
		if (s_t[ref_len][0] > 0) {
			t[ref_len][1] /= s_t[ref_len][0];
			t[ref_len][3] /= s_t[ref_len][0];
		}
			
		s_t[ref_len][1] = t[ref_len][5] + t[ref_len][6];
		if (s_t[ref_len][1] > 0) {
			t[ref_len][5] /= s_t[ref_len][1];
			t[ref_len][6] /= s_t[ref_len][1];
		}

		/* Estimate emission probabilities. */
		for (k = 0; k < ref_len; k ++) {
			s_e[k] = e[k][1] + e[k][2] + e[k][4] + e[k][8] + e[k][15];
			if (s_e[k] > 0) {
				e[k][1] /= s_e[k];
				e[k][2] /= s_e[k];
				e[k][4] /= s_e[k];
				e[k][8] /= s_e[k];
				e[k][15] /= s_e[k];
			}
			es += s_e[k];
		}
		diff = fabs(Pr - p);
		fprintf (stderr, "Pr: %g\tp: %g\tdiff: %g\ncount: %d\n", Pr, p, diff, count);
		Pr = p;
		count ++;
	}
	for (k = 0; k <= ref_len; k ++) {
		for (i = 0; i < 11; i ++) {
			fprintf (stderr, "t[%d][%d]: %g\t", k, i, t[k][i]);
		}
		fprintf (stderr, "\n");
	}

	for (k = 0; k < ref_len; k ++) {
			fprintf (stderr, "e[%d][1]: %g\te[%d][2]: %g\te[%d][4]: %g\te[%d][8]: %g\te[%d][15]: %g\n", k, e[k][1], k, e[k][2], k, e[k][4], k, e[k][8], k, e[k][15]);
	}
		
	for (i = 0; i <= ref_len; i ++) {
	    free(t[i]);
	}
	free(t);		
	for (i = 0; i < ref_len; i ++) {
		free(e[i]);
	}
	free(e);
	emission_destroy(emission, ref_len);
	transition_destroy(transition, ref_len);
} 
