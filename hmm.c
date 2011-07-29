/*
 * hmm.c: Semi-global HMM
 * Author: Mengyao Zhao
 * Create date: 2011-06-13
 * Contact: zhangmp@bc.edu
 * Last revise: 2011-07-28 
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
	
		/* Sum of the out edges of S equals to 1 */	
		matrix_array[i][9] = d/L;	/* S -> M_k+1 */
		matrix_array[i][10] = (1 - d)/(L + 1);	/* S -> I_k */
	}
	matrix_array[1][7] = 0;	/* State D_1 doesn't exist. */	

	/* k = L */
	for (i = 0; i < 11; i ++) {
		matrix_array[L][i] = 0;
	}
	matrix_array[L][1] = a*(1 - r);	/* M_k -> I_k */
	matrix_array[L][3] = r;	/* M_k -> E */
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

void emission_destroy (float** array, const int32_t L)
{
	int32_t i;
	for (i = 0; i < L; i ++) {
		free(array[i]);
	}
	free(array);
}

double forward_backward (float** transition, float** emission, char* ref, uint8_t* read, int32_t read_len, fb* f, fb* b)
{
	/* forward algorithm */
	int32_t i;	 /* iter of read */ 
	int32_t k;	 /* iter of reference */
	int32_t ref_len = strlen(ref);
	double s[read_len + 1]; /* scaling factor to avoid underflow */
	
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
		
		for (k = 2; k <= ref_len; k ++) {
			f->match[i][k] = emission[k - 1][bam1_seqi(read, i)] * (transition[k - 1][0] *
		    f->match[i - 1][k - 1] + transition[k - 1][4] * f->insertion[i - 1][k - 1] + transition[k - 1][7] *
			f->deletion[i - 1][k - 1]);
			
			f->insertion[i][k] = 0.25 * (transition[k][1] * f->match[i - 1][k] + transition[k][5] * f->insertion[i - 1][k]);
			
			f->deletion[i][k] = transition[k - 1][2] * f->match[i][k - 1] + transition[k - 1][8] * f->deletion[i][k - 1];
			
			s[i] += f->match[i][k] + f->insertion[i][k] + f->deletion[i][k];
		}
		
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
	b->match[read_len - 1][0] = 0; /* no M_0 state */
	b->insertion[read_len - 1][0] = transition[0][6];
	for (k = 1; k <= ref_len; k ++) {
		b->match[read_len - 1][k] = transition[k][3];
		b->insertion[read_len - 1][k] = transition[k][6];
	}

	/* rescale */
	for (k = 0; k <= ref_len; k ++) {
		/* b_0_E needs to be rescaled by s[read_len] */
		b->match[read_len - 1][k] = b->match[read_len - 1][k] / s[read_len - 1] / s[read_len];
		b->insertion[read_len - 1][k] = b->insertion[read_len - 1][k] / s[read_len - 1] / s[read_len];
	}
	
	for (i = read_len - 2; i > 0; i --) {
		b->match[i][ref_len] = 0.25 * transition[ref_len][1] * b->insertion[i + 1][ref_len];

		b->insertion[i][ref_len] = 0.25 * transition[ref_len][5] * b->insertion[i + 1][ref_len];

		b->deletion[i][ref_len] = 0;
	 
		for (k = ref_len  - 1; k > 0; k --) {
			b->match[i][k] = emission[k][bam1_seqi(read, i + 1)] * transition[k][0] * b->match[i + 1][k + 1] +
			0.25 * transition[k][1] * b->insertion[i + 1][k] + transition[k][2] * b->deletion[i][k + 1];

			b->insertion[i][k] = emission[k][bam1_seqi(read, i + 1)] * transition[k][4] * 
			b->match[i + 1][k + 1] + 0.25 * transition[k][5] * b->insertion[i + 1][k];

			b->deletion[i][k] = emission[k][bam1_seqi(read, i + 1)] * transition[k][7] * 
			b->match[i + 1][k + 1] + transition[k][8] * b->deletion[i][k + 1];
		}
		b->insertion[i][0] = emission[0][bam1_seqi(read, i + 1)] * transition[0][4] * b->match[i + 1][1] +
		0.25 * transition[0][5] * b->insertion[i + 1][0];
		b->match[i][0] = b->deletion[i][0] = 0; /* no M_0, D_0 states */

		/* rescale */
		for (k = 0; k <= ref_len; k ++) {
			b->match[i][k] /= s[i];
			b->insertion[i][k] /= s[i];
			b->deletion[i][k] /= s[i];
		}
	}

	b->deletion[0][ref_len] = 0;
	b->match[0][ref_len] = 0;
	b->insertion[0][ref_len] = 0.25 * transition[ref_len][5] * b->insertion[1][ref_len];
	for (k = ref_len - 1; k >= 0; k --) {
		b->match[0][k] = emission[k][bam1_seqi(read, 1)] * transition[k][0] * b->match[1][k + 1] +
		0.25 * transition[k][1] * b->insertion[i + 1][k] + transition[k][2] * b->deletion[0][k + 1];

		b->insertion[0][k] = emission[k][bam1_seqi(read, 1)] * transition[k][4] * 
		b->match[1][k + 1] + 0.25 * transition[k][5] * b->insertion[1][k];

		b->deletion[0][k] = 0;
	}

	/* rescale */
	b->final = 0;
	for (k = 1; k <= ref_len; k ++) {
		b->match[0][k] /= s[0];
		b->insertion[0][k] /= s[0];

		b->final += emission[k - 1][bam1_seqi(read, 0)] * transition[k - 1][9] * b->match[0][k] + 0.25 *
		transition[k][10] * b->insertion[0][k];
	}
	b->insertion[0][0] /= s[0];

	b->final += 0.25 * transition[0][10] * b->insertion[0][0]; 	
	{// compute the log likelihood
		double p = 1; 
		for (i = 0; i <= read_len; ++i) {
			p *= s[i];
		}
		p *= (ref_len + 1) * read_len;
		return p;
	}
}

void baum_welch (char* ref_seq, int32_t ref_len, reads* r, float df) /* 0-based coordinate */ 
{
	fprintf (stdout, "reference sequence: %s\n", ref_seq); 
	float** transition = transition_init (0.3, 0.5, 0.2, 0.5, 0.5, ref_len);
	float** emission = emission_init(ref_seq);
	double Pr = 2, diff = 1;
	int32_t i, k, j, count = 0;
	float** t = calloc (ref_len + 1, sizeof(float*)); 
	for (i = 0; i <= ref_len; i ++) {
		t[i] = calloc (11, sizeof(float));
	}		
	float** e = calloc (ref_len, sizeof(float*));
	for (i = 0; i < ref_len; i ++) {
		e[i] = calloc (16, sizeof(float));
	}

	while (diff > df && count < 100) {
		double p = 0;	/* likelihood */
		if (count > 0) {
			for (k = 0; k <= ref_len; k ++) {
				for (i = 0; i < 11; i ++) {
					transition[k][i] = t[k][i];
				}
			}
			for (k = 0; k < ref_len; k ++) {
				for (i = 0; i < 16; i ++) {
					emission[k][i] = e[k][i];
				}
			}
		}
		/* Initialize new transition and emission matrixes. */
		for (k = 0; k <= ref_len; k ++) {
			for (i = 0; i < 11; i ++) {
				t[k][i] = 0;
			}
		}
		for (k = 0; k < ref_len; k ++) {
			for (i = 0; i < 16; i ++) {
				e[k][i] = 0;
			}
		}
		float s_t[ref_len + 1][3], s_e[ref_len], s_tS = 0;
		int32_t total_hl = 0;

		/* Transition and emission matrixes training by a block of reads. */
		for (j = 0; j < r->count; j ++) {
			uint8_t* read_seq = &r->seqs[total_hl];
			total_hl += r->seq_l[j]/2 + r->seq_l[j]%2;
			int32_t read_len = r->seq_l[j];
			fprintf (stderr, "read length: %d\nread_seq: ", read_len);
			for (i = 0; i < read_len; i ++) {
				fprintf (stderr, "%d ", bam1_seqi(read_seq, i));
			}
			fprintf (stderr, "\n");

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
			p += forward_backward (transition, emission, ref_seq, read_seq, read_len, f, b);
			float temp_t[ref_len + 1][11];	/* for transition probability training */
			float temp_e[ref_len][16];	/* for emission probability training */
			for (k = 0; k < ref_len; k ++) {
				/* for transition probability training */
				for (i = 0; i < 11; i ++) {
					temp_t[k][i] = 0;
				}
				for (i = 0; i < read_len - 1; i ++) {
					temp_t[k][0] += f->match[i][k] * transition[k][0] * emission[k][bam1_seqi(read_seq, i + 1)] 
					* b->match[i + 1][k + 1];	/* M_k -> M_k+1 */

					temp_t[k][1] += 0.25 * f->match[i][k] * transition[k][1] * b->insertion[i + 1][k]; /* M_k -> I_k */

					temp_t[k][2] += f->match[i][k] * transition[k][2] * b->deletion[i][k + 1];	/* M_k -> D_k+1 */
					
					temp_t[k][4] += f->insertion[i][k] * transition[k][4] * emission[k][bam1_seqi(read_seq, i + 1)] 
					* b->match[i + 1][k + 1];	/* I_k -> M_k+1 */
			
					temp_t[k][5] += 0.25 * f->insertion[i][k] * transition[k][5] * b->insertion[i + 1][k];	/* I_k -> I_k */

					temp_t[k][7] += f->deletion[i][k] * transition[k][7] * emission[k][bam1_seqi(read_seq, i + 1)] 
					* b->match[i + 1][k + 1];	/* D_k -> M_k+1 */
		
					temp_t[k][8] += f->deletion[i][k] * transition[k][8] * b->deletion[i][k + 1];	/* D_k -> D_k+1 */
				}
				/* i = read_len - 1 */
				temp_t[k][2] += f->match[read_len - 1][k] * transition[k][2] * b->deletion[read_len - 1][k];	/* M_k -> D_k+1 */
				
				temp_t[k][3] += f->match[read_len - 1][k] * transition[k][3];	/* M_k -> E */

				temp_t[k][6] += f->insertion[read_len - 1][k] * transition[k][6];	/* I_k -> E */

				temp_t[k][8] += f->deletion[read_len - 1][k] * transition[k][8] * b->deletion[i][k + 1];	/* D_k -> D_k+1 */
				
				/* i = 0 */
				temp_t[k][9] += transition[k][9] * emission[k][bam1_seqi(read_seq, 0)] * b->match[0][k + 1];	/* S -> M_k+1 */

				temp_t[k][10] += 0.25 * transition[k][10] * b->insertion[0][k];	/* S -> I_k */
				
				/* for emission probability training */
				for (i = 0; i < 16; i ++) {
					temp_e[k][i] = 0;
				}
				for (i = 0; i < read_len; i ++) {
					temp_e[k][bam1_seqi(read_seq, i)] += f->match[i][k + 1] * b->match[i][k + 1];
					if (k == 0 && i == 0) {
						fprintf (stderr, "f->match[%d][%d]: %g\tb->match[%d][%d]: %g\n", i, k + 1, f->match[i][k + 1], i, k + 1, b->match[i][k + 1]);
					}
				}
			}

			/* for transition probability training */
			for (i = 0; i < 11; i ++) {
				temp_t[ref_len][i] = 0;
			}
			for (i = 0; i < read_len - 1; i ++) {
				 /* M_k -> I_k */
				temp_t[ref_len][1] += 0.25 * f->match[i][ref_len] * transition[ref_len][1] * b->insertion[i + 1][ref_len];
					
				/* I_k -> I_k */
				temp_t[ref_len][5] += 0.25 * f->insertion[i][ref_len] * transition[ref_len][5] * b->insertion[i + 1][ref_len];
			}
			temp_t[ref_len][3] += f->match[read_len - 1][ref_len] * transition[ref_len][3];	/* M_k -> E */
			temp_t[ref_len][6] += f->insertion[read_len - 1][ref_len] * transition[ref_len][6];	/* I_k -> E */
			temp_t[ref_len][10] += 0.25 * transition[ref_len][10] * b->insertion[0][ref_len];	/* S -> I_k */

			for (k = 0; k <= ref_len; k ++) {
				for (i = 0; i < 11; i ++) {
					t[k][i] += temp_t[k][i]/p;
				}
			}
			
			/* for emission probability training */
			for (k = 0; k < ref_len; k ++) {
				for (i = 0; i < 16; i ++) {
					e[k][i] += temp_e[k][i]/p;
				}
			}
		
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
		t[1][7] = t[1][8] = t[ref_len][0] = t[ref_len][2] = t[ref_len][4] = 
		t[ref_len][7] = t[ref_len][8] = t[ref_len][9] = 0;	
		
		/* Estimate transition probabilities. */
		s_t[0][1] = t[0][4] + t[0][5] + t[0][6];
		if (s_t[0][1] > 0) {
			t[0][4] /= s_t[0][1];
			t[0][5] /= s_t[0][1];
			t[0][6] /= s_t[0][1];
		}

		s_tS += t[0][9] + t[0][10];

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

			s_tS += t[k][9] + t[k][10];		
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

		s_tS += t[ref_len][10];
		if (s_tS > 0) {
			for (k = 0; k <= ref_len; k ++) {
				t[k][9] /= s_tS;
				t[k][10] /= s_tS;
			}
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
				fprintf (stderr, "e[%d][1]: %f\te[%d][2]: %f\te[%d][4]: %f\te[%d][8]: %f\te[%d][15]: %f\n", k, e[k][1], k, e[k][2], k, e[k][4], k, e[k][8], k, e[k][15]);
			}
		}

		p /= r->count;
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
		for (i = 0; i < 16; i ++) {
			fprintf (stderr, "*e[%d][%d]:%g*\t", k, i, e[k][i]);
		}
		fprintf (stderr, "\n");
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
