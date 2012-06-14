/*
 * hmm.c: Semi-global HMM
 * Author: Mengyao Zhao
 * Create date: 2011-06-13
 * Contact: zhangmp@bc.edu
 * Last revise: 2012-06-14 
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
   /*I_A   A     C     I_C   G     I_G   pad pad T     I_T   pad pad pad pad I_N N   
	{0.25, 1,    0,    0.25, 0,    0.25, 0,  0,  0,    0.25, 0,  0,  0,  0,  0,  0.25},	A 
	{0.25, 0,    1,    0.25, 0,    0.25, 0,  0,  0,    0.25, 0,  0,  0,  0,  0,  0.25},	C 
	{0.25, 0,    0,    0.25, 1,    0.25, 0,  0,  0,    0.25, 0,  0,  0,  0,  0,  0.25},	G 
	{0.25, 0,    0,    0.25, 0,    0.25, 0,  0,  1,    0.25, 0,  0,  0,  0,  0,  0.25},	T 
	{0.25, 0,    0,    0.25, 0.5,  0.25, 0,  0,  0.5,  0.25, 0,  0,  0,  0,  0,  0.5},	K 
	{0.25, 0.5,  0.5,  0.25, 0,    0.25, 0,  0,  0,    0.25, 0,  0,  0,  0,  0,  0.5},	M
	{0.25, 0.5,  0,    0.25, 0.5,  0.25, 0,  0,  0,    0.25, 0,  0,  0,  0,  0,  0.5},	R
	{0.25, 0,    0.5,  0.25, 0,    0.25, 0,  0,  0.5,  0.25, 0,  0,  0,  0,  0,  0.5},	Y
	{0.25, 0,    0.5,  0.25, 0.5,  0.25, 0,  0,  0,    0.25, 0,  0,  0,  0,  0,  0.5},	S
	{0.25, 0,    0.33, 0.25, 0.33, 0.25, 0,  0,  0.33, 0.25, 0,  0,  0,  0,  0,  0.75},	B
	{0.25, 0.33, 0.33, 0.25, 0.33, 0.25, 0,  0,  0,    0.25, 0,  0,  0,  0,  0,  0.75},	V
	{0.25, 0.33, 0.33, 0.25, 0,    0.25, 0,  0,  0.33, 0.25, 0,  0,  0,  0,  0,  0.75},	H
	{0.25, 0.33, 0,    0.25, 0.33, 0.25, 0,  0,  0.33, 0.25, 0,  0,  0,  0,  0,  0.75},	D
	{0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0,  0,  0.25, 0.25, 0,  0,  0,  0,  0,  1},	N
	{0.25, 0,    0,    0.25, 0,    0.25, 0,  0,  0,    0.25, 0,  0,  0,  0,  0,  0},	X
	{0.25, 0,    0,    0.25, 0,    0.25, 0,  0,  0,    0.25, 0,  0,  0,  0,  0,  0},	pad */
	int32_t i;
	int32_t ref_len = strlen(ref);
	double ** array = (double**)calloc(ref_len + 1, sizeof(double*));
	
	array[0] = (double*)calloc(16, sizeof(double));
	array[0][0] = array[0][3] = array[0][5] = array[0][9] = 0.25;
	array[0][14] = 3e-8;
	
	for (i = 0; i < ref_len; i ++) {
		array[i + 1] = (double*)calloc(16, sizeof(double));
		array[i + 1][0] = array[i + 1][3] = array[i + 1][5] = array[i + 1][9] = 0.25;
		array[i + 1][14] = 3e-8;
		switch (ref[i]) {
			case 'A':
			case 'a':
				array[i + 1][1] = 1; array[i + 1][15] = 0.25; array[i + 1][2] = array[i + 1][4] = array[i + 1][8] = 0.001;	
				break;
			case 'C':
			case 'c':
				array[i + 1][2] = 1; array[i + 1][15] = 0.25; array[i + 1][1] = array[i + 1][4] = array[i + 1][8] = 0.001;	
				break;
			case 'G':
			case 'g':
				array[i + 1][4] = 1; array[i + 1][15] = 0.25; array[i + 1][1] = array[i + 1][2] = array[i + 1][8] = 0.001;	
				break;
			case 'T':
			case 't':
				array[i + 1][8] = 1; array[i + 1][15] = 0.25; array[i + 1][1] = array[i + 1][2] = array[i + 1][4] = 0.001;	
				break;
			case 'N': /* any */
			case 'n':
				array[i + 1][1] = array[i + 1][2] = array[i + 1][4] = array[i + 1][8] = 0.25; array[i + 1][15] = 1;	
				break;
			case 'X': /* any, x mask on Y chromosome */
			case 'x':
				array[i + 1][1] = array[i + 1][2] = array[i + 1][4] = array[i + 1][8] = array[i + 1][15] = 0.001;	
				break;
			case 'K': /* G or T */
			case 'k':
				array[i + 1][4] = array[i + 1][8] = array[i + 1][15] = 0.5; array[i + 1][1] = array[i + 1][2] = 0.001;	
				break;
			case 'M': /* A or C */
			case 'm':
				array[i + 1][1] = array[i + 1][2] = array[i + 1][15] = 0.5; array[i + 1][4] = array[i + 1][8] = 0.001;	
				break;
			case 'R': /* A or G */
			case 'r':
				array[i + 1][1] = array[i + 1][4] = array[i + 1][15] = 0.5; array[i + 1][2] = array[i + 1][8] = 0.001;	
				break;
			case 'Y': /* C or T */
			case 'y':
				array[i + 1][2] = array[i + 1][8] = array[i + 1][15] = 0.5; array[i + 1][1] = array[i + 1][4] = 0.001;	
				break;
			case 'S': /* A or T */
			case 's':
				array[i + 1][2] = array[i + 1][4] = array[i + 1][15] = 0.5; array[i + 1][1] = array[i + 1][8] = 0.001;	
				break;
			case 'B': /* C or G or T */
			case 'b':
				array[i + 1][2] = array[i + 1][4] = array[i + 1][8] = 0.33; array[i + 1][15] = 0.75; array[i + 1][1] = 0.001;	
				break;
			case 'V': /* A or C or G */
			case 'v':
				array[i + 1][1] = array[i + 1][2] = array[i + 1][4] = 0.33; array[i + 1][15] = 0.75; array[i + 1][8] = 0.001;	
				break;
			case 'H': /* A or C or T */
			case 'h':
				array[i + 1][1] = array[i + 1][2] = array[i + 1][8] = 0.33; array[i + 1][15] = 0.75; array[i + 1][4] = 0.001;	
				break;
			case 'D': /* A or G or T */
			case 'd':
				array[i + 1][1] = array[i + 1][4] = array[i + 1][8] = 0.33; array[i + 1][15] = 0.75; array[i + 1][2] = 0.001;	
				break;
			default:
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
	for (i = 0; i <= L; i ++) {
		free(array[i]);
	}
	free(array);
}

double forward_backward (double** transition, 
						 double** emission, 
						 int32_t ref_begin,
						 int32_t ref_len,	// region size 
						 uint8_t* read, 
						 int32_t read_len, 
						 double** f, 
						 double** b, 
						 double* s, 
						 int32_t bw) {
	
	#define set_u(u, b, i, k) { int32_t x=(i)-(b); x=x>0?x:0; (u)=((k)-x+1)*3; }
 
	/*-------------------*
	 * forward algorithm *
	 *-------------------*/
	int32_t i;	 /* iter of read */ 
	int32_t k;	 /* iter of reference */
	int32_t temp1, temp, x, u; 
	int32_t beg = ref_begin - bw > 0 ? ref_begin - bw : 0, end = ref_len < ref_begin + bw ? ref_len : ref_begin + bw;
	double f_final, b_final;
//	 double pp = 0;	// Debug: posterior probability of each state 
// f[0]
	temp1 = bam1_seqi(read, 0);
	temp = temp1 + pow(-1, temp1%2);
	set_u(u, bw, 0, beg);
	f[0][u + 1] = transition[beg][10] * emission[beg][temp];	// 1: insertion
	s[0] = f[0][u + 1]; 	// 1: insertion

	for (k = beg + 1; k <= end; k ++) {
		set_u(u, bw, 0, k);
		f[0][u] = emission[k][bam1_seqi(read, 0)] * transition[k - 1][9];	// 0: match
		f[0][u + 1] = transition[k][10] * emission[k][temp];	// 1: insertion
		s[0] += f[0][u] + f[0][u + 1];
	}

	/* rescale */
	for (k = beg; k <= end; k ++) {
		set_u(u, bw, 0, k);
		f[0][u] /= s[0];	// 0: match
		f[0][u + 1] /= s[0];	// 1: insertion
	}

//	f[i]
	for (i = 1; i < read_len; i ++) {
		int32_t v, w;
		temp1 = bam1_seqi(read, i);
		temp = temp1 + pow(-1, temp1%2);

		beg = ref_begin + i - bw > 0 ? ref_begin + i - bw : 0;
		set_u(u, bw, i, beg);
		set_u(v, bw, i - 1, beg);
		set_u(w, bw, i - 1, beg - 1);
		f[i][u + 1] = transition[beg][5] * emission[beg][temp] * f[i - 1][v + 1]; /* f_i_I0; i = 2 ~ l */
		
		set_u(w, bw, i, beg + 1);
		f[i][w] = emission[beg + 1][bam1_seqi(read, i)] * transition[beg][4] * f[i - 1][v + 1]; /* f_i_M1; i = 1 ~ l */
		
		set_u(v, bw, i - 1, beg + 1);
		f[i][w + 1] = emission[beg + 1][temp] * (transition[beg + 1][1] * f[i - 1][v] + transition[beg + 1][5] * f[i - 1][v + 1]); /* f_i_I1; i = 2 ~ l */

		s[i] = f[i][u + 1] + f[i][w] + f[i][w + 1];

		end = ref_len;	
		x = ref_begin + i + bw; end = end < x ? end : x; //	band end
		for (k = beg + 2; k < end; k ++) {
			set_u(u, bw, i, k);
			set_u(v, bw, i - 1, k - 1);
			f[i][u] = emission[k][bam1_seqi(read, i)] * (transition[k - 1][0] *	// 0: match
		    f[i - 1][v] + transition[k - 1][4] * f[i - 1][v + 1] + transition[k - 1][7] * f[i - 1][v + 2]);
			
			set_u(w, bw, i - 1, k);
			f[i][u + 1] = emission[k][temp] * (transition[k][1] * f[i - 1][w] + transition[k][5] * f[i - 1][w + 1]);	// 1: insertion
			set_u(v, bw, i, k - 1);
			f[i][u + 2] = transition[k - 1][2] * f[i][v] + transition[k - 1][8] * f[i][v + 2];	// 2: deletion
			
			s[i] += f[i][u] + f[i][u + 1] + f[i][u + 2];
		}

		set_u(u, bw, i, end);
		set_u(v, bw, i - 1, end - 1);
		f[i][u] = emission[end][bam1_seqi(read, i)] * (transition[end - 1][0] *	// 0: match
		f[i - 1][v] + transition[end - 1][4] * f[i - 1][v + 1] + transition[end - 1][7] * f[i - 1][v + 2]);
		
		set_u(w, bw, i - 1, end);
		f[i][u + 1] = emission[end][temp] * (transition[end][1] * f[i - 1][w] + transition[end][5] * f[i - 1][w + 1]);	// 1: insertion
		
		s[i] += f[i][u] + f[i][u + 1];
		
		/* rescale */
		for (k = beg; k <= end; k ++) {
			set_u(u, bw, i, k);
			f[i][u] /= s[i];	// 0: match
			f[i][u + 1] /= s[i];	// 1: insertion
			f[i][u + 2] /= s[i];	// 2: deletion
		}
	}

	/* sum of all forward path */
	f_final = 0;
	beg = ref_begin + read_len - bw; end = ref_begin + read_len + bw < ref_len ? ref_begin + read_len + bw : ref_len;
	for (k = beg; k <= end; k ++) {
		set_u(u, bw, read_len - 1, k)
		f_final += transition[k][3] * f[read_len - 1][u] + transition[k][6] * f[read_len - 1][u + 1];
	}
	
	s[read_len] = f_final;
//	f_final /= s[read_len];	// Debug

	/*--------------------*
	 * backword algorithm *
	 *--------------------*/
	set_u(u, bw, read_len - 1, beg);
	b[read_len - 1][u + 1] = transition[beg][6] / s[read_len];	// 1: insertion
	for (k = beg; k <= end; k ++) {
		set_u(u, bw, read_len - 1, k);
		b[read_len - 1][u] = transition[k][3] / s[read_len];	// 0: match
		b[read_len - 1][u + 1] = transition[k][6] / s[read_len];	// 1: insertion
	}

	/* rescale */
	if (read_len > 1) {	
		for (k = beg; k <= end; k ++) {
			/* b_0_E needs to be rescaled by s[read_len] */
			set_u(u, bw, read_len - 1, k);
			b[read_len - 1][u] /= s[read_len - 1];	// 0: match
			b[read_len - 1][u + 1] /= s[read_len - 1];	// 1: insertion
			/* Debug: posterior probability */
//			 pp += b[read_len - 1][3*k] * f[read_len - 1][3*k] + b[read_len - 1][3*k + 1] * f[read_len - 1][3*k + 1];	// Debug
		}
//		pp *= s[read_len - 1];	// Debug
//		fprintf (stderr, "pp: %f\n", pp);	// Debug
	
		for (i = read_len - 2; i > 0; i --) {
			int32_t v, w, y;
			temp1 = bam1_seqi(read, i + 1);
			temp = temp1 + pow(-1, temp1%2);
			end = ref_len;
			x = ref_begin + i + bw; end = end < x ? end : x; //	band end

			set_u(u, bw, i, end);
			set_u(v, bw, i + 1, end);
			b[i][u] = transition[end][1] * emission[end][temp] * b[i + 1][v + 1];	// 0: match

			b[i][u + 1] = transition[end][5] * emission[end][temp] * b[i + 1][v + 1];	// 1: insertion
			
			set_u(u, bw, i, end - 1);
			set_u(w, bw, i + 1, end - 1);
			b[i][u] = emission[end][bam1_seqi(read, i + 1)] * transition[end - 1][0] * b[i + 1][v] +
			emission[end - 1][temp] * transition[end - 1][1] * b[i + 1][w + 1];	/* 0: match; no D_L state */

			b[i][u + 1] = emission[end][bam1_seqi(read, i + 1)] * transition[end - 1][4] * 
			b[i + 1][v] + transition[end - 1][5] * emission[end - 1][temp] * b[i + 1][w + 1];	// 1: insertion

			b[i][u + 2] = emission[end][bam1_seqi(read, i + 1)] * transition[end - 1][7] * b[i + 1][v];	/* 2: deletion; no D_L state */
			
			beg = 0; 
			x = ref_begin + i - bw; beg = beg > x ? beg : x;
			for (k = end - 2; k > beg; k --) {
				set_u(u, bw, i, k);
				set_u(v, bw, i + 1, k + 1);
				set_u(w, bw, i + 1, k);
				set_u(y, bw, i, k + 1);
				b[i][u] = emission[k + 1][bam1_seqi(read, i + 1)] * transition[k][0] * b[i + 1][v] +
				transition[k][1] * emission[k][temp] * b[i + 1][w + 1] + transition[k][2] * b[i][y + 2];	// 0: match

				b[i][u + 1] = emission[k + 1][bam1_seqi(read, i + 1)] * transition[k][4] * 
				b[i + 1][v] + transition[k][5] * emission[k][temp] * b[i + 1][w + 1];	// 1: insertion

				b[i][u + 2] = emission[k + 1][bam1_seqi(read, i + 1)] * transition[k][7] * 
				b[i + 1][v] + transition[k][8] * b[i][y + 2];	// 2: deletion
			}

			set_u(u, bw, i, beg);
			set_u(v, bw, i + 1, beg + 1);
			set_u(w, bw, i + 1, beg);
			b[i][u + 1] = emission[beg + 1][bam1_seqi(read, i + 1)] * transition[beg][4] * b[i + 1][v] +
			transition[beg][5] * emission[beg][temp] * b[i + 1][w + 1];	// 1: insertion

			/* rescale */
//			pp = 0;	// Debug 
			for (k = beg; k <= end; k ++) {
				set_u(u, bw, i, k);
				b[i][u] /= s[i];	// 0: match
				b[i][u + 1] /= s[i];	// 1: insertion
				b[i][u + 2] /= s[i];	// 2: deletion
//				pp += b[i][3*k] * f[i][3*k] + b[i][3*k + 1] * f[i][3*k + 1]; // Debug
			}
//			pp *= s[i];	// Debug
//			fprintf (stderr, "pp: %f\n", pp);	// Debug
		}

		temp1 = bam1_seqi(read, 1);
		temp = temp1 + pow(-1, temp1%2);
		end = ref_len < ref_begin + bw ? ref_len : ref_begin + bw;
		set_u(u, bw, 0, end);
		set_u(v, bw, 1, end);
		b[0][u] = transition[end][1] * emission[end][temp] * b[1][v + 1];	// 0: match
		b[0][u + 1] = transition[end][5] * emission[end][temp] * b[1][v + 1];	// 1: insertion

		set_u(w, bw, 0, end - 1);
		set_u(y, bw, 1, end - 1);
		b[0][w] = emission[end][bam1_seqi(read, 1)] * transition[end - 1][0] * b[1][v] +
		transition[end - 1][1] * emission[end - 1][temp] * b[1][y + 1];	/* 0: match; no D_L state */

		b[0][w + 1] = emission[end][bam1_seqi(read, 1)] * transition[end - 1][4] * 
		b[1][v] + transition[end - 1][5] * emission[end - 1][temp] * b[1][y + 1];	// 1: insertion

		beg = ref_begin - bw > 0 ? ref_begin - bw : 0;
		for (k = end - 2; k >= beg; k --) {
			set_u(u, bw, 0, k);
			set_u(v, bw, 1, k + 1);
			set_u(w, bw, 1, k);
			set_u(y, bw, 0, k + 1);
			b[0][u] = emission[k + 1][bam1_seqi(read, 1)] * transition[k][0] * b[1][v] +
			transition[k][1] * emission[k][temp] * b[1][w + 1] + transition[k][2] * b[0][y + 2];	// 0: match

			b[0][u + 1] = emission[k + 1][bam1_seqi(read, 1)] * transition[k][4] * 
			b[1][v] + transition[k][5] * emission[k][temp] * b[1][w + 1];	// 1: insertion
		}
	}
	/* rescale */
//	pp = 0; // Debug
	b_final = 0;
	temp1 = bam1_seqi(read, 0);
	temp = temp1 + pow(-1, temp1%2);

	beg = ref_begin - bw > 0 ? ref_begin - bw : 0;
	end = ref_len < ref_begin + bw ? ref_len : ref_begin + bw;
	for (k = beg + 1; k <= end; k ++) {
		set_u(u, bw, 0, k);
		b[0][u] /= s[0];	// 0: match
		b[0][u + 1] /= s[0];	// 1: insertion
//		pp += b[0][3*k] * f[0][3*k] + b[0][3*k + 1] * f[0][3*k + 1];	// Debug: 0: match

		b_final += emission[k][bam1_seqi(read, 0)] * transition[k - 1][9] * b[0][u] + 
		transition[k][10] * emission[k][temp] * b[0][u + 1];
	}
	set_u(u, bw, 0, 0);
	b[0][u + 1] /= s[0];	// 1: insertion
//	pp += b[0][1] * f[0][1];	// Debug: 1: insertion
//	pp *= s[0]; // Debug
	b_final += transition[beg][10] * emission[beg][temp] * b[0][u + 1];
//	fprintf (stderr, "pp: %f\n", pp);	// Debug
//	fprintf (stderr, "b_final: %g\n", b_final);	// Debug: b_final should equal to 1 

	// Debug: posterior probability for transition (each edge) 
/*	{
		double pp_t = 0;
		for (i = 0; i < read_len - 1; i ++) {
			pp_t = 0;
			temp1 = bam1_seqi(read, i + 1);
			temp = temp1 + pow(-1, temp1%2);
			for (k = 2; k < ref_len; k ++) {
				pp_t += f[i][3*k + 2] * transition[k][7] * emission[k + 1][bam1_seqi(read, i + 1)] * b[i + 1][3*(k + 1)];
			}
			for (k = 1; k < ref_len; k ++) {
				pp_t += f[i][3*k] * transition[k][0] * emission[k + 1][bam1_seqi(read, i + 1)] * b[i + 1][3*(k + 1)];
				pp_t += f[i][3*k + 1] * transition[k][4] * emission[k + 1][bam1_seqi(read, i + 1)] * b[i + 1][3*(k + 1)];
				pp_t += emission[k][temp] * f[i][3*k] * transition[k][1] * b[i + 1][3*k + 1];
				pp_t += emission[k][temp] * f[i][3*k + 1] * transition[k][5] * b[i + 1][3*k + 1];
			}
			pp_t += f[i][1] * transition[0][4] * emission[1][bam1_seqi(read, i + 1)] * b[i + 1][3];
			pp_t += emission[ref_len][temp] * f[i][3*ref_len] * transition[ref_len][1] * b[i + 1][3*ref_len + 1]; 
			pp_t += emission[0][temp] * f[i][1] * transition[0][5] * b[i + 1][1];
			pp_t += emission[ref_len][temp] * f[i][3*ref_len + 1] * transition[ref_len][5] * b[i + 1][3*ref_len + 1];
			fprintf (stderr, "pp_t: %g\n", pp_t);
		}
	}*/ // Debug
		
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

void baum_welch (double** transition, double** emission, char* ref_seq, int32_t ref_len, reads* r, double df) /* 0-based coordinate */ 
{
	double Pr = 10e100, diff = 1;
	int32_t i, k, j, count = 0;
	double** t = calloc (ref_len + 1, sizeof(double*));
	double** e = calloc (ref_len + 1, sizeof(double*));
	for (i = 0; i <= ref_len; i ++) {
		t[i] = calloc (16, sizeof(double));
		e[i] = calloc (16, sizeof(double));
	}		
	while (diff > df && count < 10) {
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
			for (i = 0; i < 16; i ++) {
				e[k][i] = 0;
			}
		}
		double** s_t = calloc (ref_len + 1, sizeof(double*));
		for (k = 0; k < ref_len + 1; k ++) {
			s_t[k] = calloc (3, sizeof(double));
		}
		double** s_e = calloc (ref_len + 1, sizeof(double*));
		for (k = 0; k < ref_len + 1; k ++) {
			s_e[k] = calloc (2, sizeof(double));
		}
		int32_t total_hl = 0;

		/* Transition and emission matrixes training by a block of reads. */
		for (j = 0; j < r->count; j ++) {
			uint8_t* read_seq = &r->seqs[total_hl];
			total_hl += r->seq_l[j]/2 + r->seq_l[j]%2;
			int32_t read_len = r->seq_l[j];

			int32_t bw = 100, bw2;
		//	if (bw < abs(ref_len - read_len) bw = abs(ref_len - read_len);
			bw2 = bw * 2 + 1;

			int32_t temp1;
			double temp;

			double** f = (double**)calloc(read_len, sizeof(double*));
			double** b = (double**)calloc(read_len, sizeof(double*));
			for (i = 0; i < read_len; ++i) {
			//	f[i] = (double*)calloc(3*(ref_len + 1), sizeof(double));
			//	b[i] = (double*)calloc(3*(ref_len + 1), sizeof(double));
				f[i] = (double*)calloc(3*bw2, sizeof(double));
				b[i] = (double*)calloc(3*bw2, sizeof(double));
			}
			double* s = (double*)calloc(read_len + 1, sizeof(double));			

			p += forward_backward (transition, emission, ref_begin, read_seq, read_len, f, b, s, bw);

			for (k = 0; k < ref_len; k ++) {
				for (i = 0; i < read_len - 1; i ++) {
					temp1 = bam1_seqi(read_seq, i + 1);
					temp = emission[k][temp1 + (int32_t)pow(-1, temp1%2)];
					t[k][0] += f[i][3*k] * transition[k][0] * emission[k + 1][bam1_seqi(read_seq, i + 1)] 
					* b[i + 1][3*(k + 1)];	/* M_k -> M_k+1 */

					t[k][1] += f[i][3*k] * transition[k][1] * temp * b[i + 1][3*k + 1]; /* M_k -> I_k */

					t[k][2] += f[i][3*k] * transition[k][2] * b[i][3*(k + 1) + 2] * s[i];	/* M_k -> D_k+1 */
					
					t[k][4] += f[i][3*k + 1] * transition[k][4] * emission[k + 1][bam1_seqi(read_seq, i + 1)] 
					* b[i + 1][3*(k + 1)];	/* I_k -> M_k+1 */
			
					t[k][5] += f[i][3*k + 1] * transition[k][5] * temp * b[i + 1][3*k + 1];	/* I_k -> I_k */

					t[k][7] += f[i][3*k + 2] * transition[k][7] * emission[k + 1][bam1_seqi(read_seq, i + 1)] 
					* b[i + 1][3*(k + 1)];	/* D_k -> M_k+1 */
		
					t[k][8] += f[i][3*k + 2] * transition[k][8] * b[i][3*(k + 1) + 2] * s[i];	/* D_k -> D_k+1 */
				}

				/* i = read_len - 1 */
				t[k][2] += f[read_len - 1][3*k] * transition[k][2] * b[read_len - 1][3*(k + 1) + 2] * s[read_len - 1];	/* M_k -> D_k+1 */
				
				t[k][8] += f[read_len - 1][3*k + 2] * transition[k][8] * b[read_len - 1][3*(k + 1) + 2] * s[read_len - 1];	/* D_k -> D_k+1 */
				
				for (i = 0; i < read_len; i ++) {
					e[k + 1][bam1_seqi(read_seq, i)] += f[i][3*(k + 1)] * b[i][3*(k + 1)] * s[i];	/* M_k+1 */
					
					temp1 = bam1_seqi(read_seq, i);
					e[k][temp1 + (int32_t)pow(-1, temp1%2)] += f[i][3*k + 1] * b[i][3*k + 1] * s[i];	/* I_k */
				}
			}
			for (i = 0; i < read_len - 1; i ++) {
				temp1 = bam1_seqi(read_seq, i + 1);
				temp = emission[ref_len][temp1 + (int32_t)pow(-1, temp1%2)];
				 /* M_k -> I_k */
				t[ref_len][1] += f[i][3*ref_len] * transition[ref_len][1] * temp * b[i + 1][3*ref_len + 1];
					
				/* I_k -> I_k */
				t[ref_len][5] += f[i][3*ref_len + 1] * transition[ref_len][5] * temp * b[i + 1][3*ref_len + 1];
			}
			for (i = 0; i < read_len; i ++) {
				temp1 = bam1_seqi(read_seq, i);
				e[ref_len][temp1 + (int32_t)pow(-1, temp1%2)] += f[i][3*ref_len + 1] * b[i][3*ref_len + 1] * s[i];	/* I_k */
			}
			free(s);
			for (i = 0; i < read_len; ++i) {
				free(f[i]);
				free(b[i]);
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
		/* Match states */
		s_e[0][1] = e[0][0] + e[0][3] + e[0][5] + e[0][9] + e[0][14];
		if (s_e[0][1] > 0) {
			e[0][0] /= s_e[0][1];
			e[0][3] /= s_e[0][1];
			e[0][5] /= s_e[0][1];
			e[0][9] /= s_e[0][1];
			e[0][14] /= s_e[0][1];
		}
		for (k = 1; k <= ref_len; k ++) {
			s_e[k][0] = e[k][1] + e[k][2] + e[k][4] + e[k][8] + e[k][15];
			if (s_e[k][0] > 0) {
				e[k][1] /= s_e[k][0];
				e[k][2] /= s_e[k][0];
				e[k][4] /= s_e[k][0];
				e[k][8] /= s_e[k][0];
				e[k][15] /= s_e[k][0];
			}
			s_e[k][1] = e[k][0] + e[k][3] + e[k][5] + e[k][9] + e[k][14];
			if (s_e[k][1] > 0) {
				e[k][0] /= s_e[k][1];
				e[k][3] /= s_e[k][1];
				e[k][5] /= s_e[k][1];
				e[k][9] /= s_e[k][1];
				e[k][14] /= s_e[k][1];
			}
		}

		for (k = 0; k < ref_len + 1; k ++) free(s_e[k]);
		free(s_e);
		for (k = 0; k < ref_len + 1; k ++) free(s_t[k]);
		free(s_t);

		diff = fabs(Pr - p);
		Pr = p;
		count ++;
	}
	for (k = 0; k <= ref_len; k ++) {
		for (i = 0; i < 16; i ++) {
			transition[k][i] = t[k][i];
	//		fprintf(stderr, "t[%d][%d]: %g\t", k, i, t[k][i]);
			emission[k][i] = e[k][i];
		}
	//	fprintf(stderr, "\n");
	}
	for (i = 0; i <= ref_len; i ++) {
	    free(t[i]);
		free(e[i]);
	}
	free(t);		
	free(e);
} 
