/*
 * hmm.c: Semi-global HMM
 * Author: Mengyao Zhao
 * Create date: 2011-06-13
 * Contact: zhangmp@bc.edu
 * Last revise: 2014-04-08 
 */

#include <math.h>
#include "bam.h"
#include "hmm.h"

#define set_u(u, b, i, k) (u)=((k)-(i)+(b))*3;

double** transition_init (const double a, const double b, const double r, const double c, const double d, const int32_t L)
{
	double** matrix_array = calloc (L + 1, sizeof(double*));
	int32_t i;
	for (i = 0; i < L + 1; i ++) {
		matrix_array[i] = calloc (16, sizeof(double));
	}

	/* k = 0: inseart before the reference */	
	matrix_array[0][4] = 1 - r - c*r;	/* I_k -> M_k+1 */
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

		/* Sum of the followin:g 3 lines equals to 1. */
		matrix_array[i][4] = 1 - r - c*r;	/* I_k -> M_k+1 */
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
	matrix_array[L - 1][4] = 1 - r - c*r;	/* I_k -> M_k+1 */
	matrix_array[L - 1][5] = c*r;	/* I_k -> I_k */ 
	matrix_array[L - 1][6] = r;	/* I_k -> E */

	/* Sum of the following 2 lines equals to 1. */
	matrix_array[L - 1][7] = 1;	/* D_k -> M_k+1 */
	matrix_array[L - 1][8] = 0;		/* D_k -> D_k+1, no D_L state */

	/* Sum of the out edges of S equals to 1 */	
	matrix_array[L - 1][9] = d/L;	/* S -> M_k+1 */
	matrix_array[L - 1][10] = (1 - d)/(L + 1);	/* S -> I_k */
	
	/* k = L */
	matrix_array[L][1] = 1 - r;	/* M_k -> I_k */
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

double** emission_init (char* ref, 
						const double x, 	// <= 0.25
						const double y, 	// <= 1
						const double z) 	// <= 0.33
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
	array[0][0] = array[0][3] = array[0][5] = array[0][9] = x;
	array[0][14] = 1 - 4*x;
	
	for (i = 0; i < ref_len; i ++) {
		array[i + 1] = (double*)calloc(16, sizeof(double));
		array[i + 1][0] = array[i + 1][3] = array[i + 1][5] = array[i + 1][9] = x;
		array[i + 1][14] = 1 - 4*x;
		switch (ref[i]) {
			case 'A':
			case 'a':
				array[i + 1][1] = y; array[i + 1][15] = array[i + 1][2] = array[i + 1][4] = array[i + 1][8] = (1 - y)/4;	
				break;
			case 'C':
			case 'c':
				array[i + 1][2] = y; array[i + 1][15] = array[i + 1][1] = array[i + 1][4] = array[i + 1][8] = (1 - y)/4;	
				break;
			case 'G':
			case 'g':
				array[i + 1][4] = y; array[i + 1][15] = array[i + 1][1] = array[i + 1][2] = array[i + 1][8] = (1 - y)/4;	
				break;
			case 'T':
			case 't':
				array[i + 1][8] = y; array[i + 1][15] = array[i + 1][1] = array[i + 1][2] = array[i + 1][4] = (1 - y)/4;	
				break;
			case 'N': /* any */
			case 'n':
				array[i + 1][1] = array[i + 1][2] = array[i + 1][4] = array[i + 1][8] = array[i + 1][15] = 0.2;	
				break;
			case 'X': /* any, x mask on Y chromosome */
			case 'x':
				array[i + 1][1] = array[i + 1][2] = array[i + 1][4] = array[i + 1][8] = array[i + 1][15] = 0.2;	
				break;
			case 'K': /* G or T */
			case 'k':
				array[i + 1][4] = array[i + 1][8] = array[i + 1][15] = z; array[i + 1][1] = array[i + 1][2] = (1 - 3*z)/2;	
				break;
			case 'M': /* A or C */
			case 'm':
				array[i + 1][1] = array[i + 1][2] = array[i + 1][15] = z; array[i + 1][4] = array[i + 1][8] = (1 - 3*z)/2;	
				break;
			case 'R': /* A or G */
			case 'r':
				array[i + 1][1] = array[i + 1][4] = array[i + 1][15] = z; array[i + 1][2] = array[i + 1][8] = (1 - 3*z)/2;	
				break;
			case 'Y': /* C or T */
			case 'y':
				array[i + 1][2] = array[i + 1][8] = array[i + 1][15] = z; array[i + 1][1] = array[i + 1][4] = (1 - 3*z)/2;	
				break;
			case 'S': /* A or T */
			case 's':
				array[i + 1][2] = array[i + 1][4] = array[i + 1][15] = z; array[i + 1][1] = array[i + 1][8] = (1 - 3*z)/2;	
				break;
			case 'B': /* C or G or T */
			case 'b':
				array[i + 1][2] = array[i + 1][4] = array[i + 1][8] = array[i + 1][15] = x; array[i + 1][1] = 1 - 4*x;	
				break;
			case 'V': /* A or C or G */
			case 'v':
				array[i + 1][1] = array[i + 1][2] = array[i + 1][4] = array[i + 1][15] = x; array[i + 1][8] = 1 - 4*x;	
				break;
			case 'H': /* A or C or T */
			case 'h':
				array[i + 1][1] = array[i + 1][2] = array[i + 1][8] = array[i + 1][15] = x; array[i + 1][4] = 1 - 4*x;	
				break;
			case 'D': /* A or G or T */
			case 'd':
				array[i + 1][1] = array[i + 1][4] = array[i + 1][8] = array[i + 1][15] = x; array[i + 1][2] = 1 - 4*x;	
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
						 int32_t ref_begin,	// relative read mapping location in window 1_based
						 int32_t window_len,	// window size 
						 uint8_t* read, 
						 int32_t read_len, 
						 double** f, 
						 double** b, 
						 double* s, 
						 int32_t bw) {

	/*-------------------*
	 * forward algorithm *
	 *-------------------*/
	int32_t i;	 /* iter of read */ 
	int32_t k;	 /* iter of reference */
	int32_t temp1, temp, x, u, v, w, y, bw2 = 3*(2*bw + 1), even; 
	int32_t beg = ref_begin - bw > 0 ? ref_begin - bw : 0, end = window_len < ref_begin + bw ? window_len : ref_begin + bw;
	double f_final, b_final;

	for (i = 0; i < read_len; ++i) fprintf(stderr, "read[%d]: %d\t", i, bam1_seqi(read, i));
	fprintf(stderr, "\nbeg: %d\n", beg);

/*	for (i = 0; i <= window_len; ++i) {
		for (k = 0; k < 16; ++k) fprintf(stderr, "e[%d][%d]: %g\t", i, k, emission[i][k]);
		fprintf(stderr, "\n");
	}*/

#ifdef VERBOSE_DEBUG
	double pp = 0;	// Debug: posterior probability of each state 
#endif
	
// f[0]
	temp1 = bam1_seqi(read, 0);
	temp = temp1 + pow(-1, temp1%2);
	set_u(u, bw, 0, 0 - ref_begin);
	s[0] = 0;

	if (beg == 0) {
		f[0][u + 1] = transition[0][10] * emission[0][temp];	// 1: insertion
//fprintf(stderr, "f[0][%d]: %g\tt[0][10]: %g\te[0][%d]: %g\n", u + 1, f[0][u + 1], transition[0][10], temp, emission[0][temp]);
		s[0] += f[0][u + 1]; 	// 1: insertion
	}

	if (beg <= 1) {
		set_u(u, bw, 0, 1 - ref_begin);
		f[0][u] = transition[0][9] * emission[1][temp1];	// 0: match
		f[0][u + 1] = transition[1][10] * emission[1][temp];	// 1: insertion
		s[0] += f[0][u] + f[0][u + 1];
	}

	even = beg > 2 ? beg : 2;
	for (k = even; k <= end; k ++) {
		set_u(u, bw, 0, k - ref_begin);
		set_u(v, bw, 0, k - 1 - ref_begin);
		f[0][u] = transition[k - 1][9] * emission[k][temp1];	// 0: match
		f[0][u + 1] = transition[k][10] * emission[k][temp];	// 1: insertion
		if (k < window_len) f[0][u + 2] = transition[k - 1][2] * f[0][v] + transition[k - 1][8] * f[0][v + 2];	//	2: deletion
		s[0] += f[0][u] + f[0][u + 1] + f[0][u + 2];
	}
	
	/* rescale */
	for (k = beg; k <= end; k ++) {
		set_u(u, bw, 0, k - ref_begin);
//fprintf(stderr, "s[0]: %g\n", s[0]);
		f[0][u] /= s[0];	// 0: match
		f[0][u + 1] /= s[0];	// 1: insertion
		f[0][u + 2] /= s[0];	// 2: deletion
	}

//	f[i]
	for (i = 1; i < read_len; i ++) {
		int32_t even;		
		temp1 = bam1_seqi(read, i);
		temp = temp1 + pow(-1, temp1%2);
		s[i] = 0;

		x = ref_begin + i - bw; beg = x > 0 ? x : 0;
		if (beg == 0) {
			set_u(u, bw, i, 0 - ref_begin);
			set_u(v, bw, i - 1, 0 - ref_begin);
			f[i][u + 1] = f[i - 1][v + 1] * transition[0][5] * emission[0][temp];	// f_i_I0	
			s[i] += f[i][u + 1];
//fprintf(stderr, "beg0\ts[%d]: %g\tf[%d][%d]: %g\n", i, s[i], i - 1, v + 1, f[i - 1][v + 1]);
		}
		
		if (beg <= 1) {	
			set_u(v, bw, i - 1, 0 - ref_begin);
			set_u(w, bw, i, 1 - ref_begin);
			f[i][w] =  transition[0][4] * emission[1][temp1] * f[i - 1][v + 1]; // f_i_M1; i = 1 ~ l 
		
			set_u(v, bw, i - 1, 1 - ref_begin);
			f[i][w + 1] = emission[1][temp] 
			* (transition[1][1] * f[i - 1][v] + transition[1][5] * f[i - 1][v + 1]); /* f_i_I1; i = 2 ~ l */
			s[i] += f[i][w] + f[i][w + 1];
		}
		
		x = ref_begin + i + bw; end = window_len < x ? window_len : x; //	band end
		even = beg > 2 ? beg : 2;
		for (k = even; k <= end; k ++) {
			set_u(u, bw, i, k - ref_begin);
			set_u(v, bw, i - 1, k - 1 - ref_begin);
			if (v >= 0) f[i][u] = emission[k][temp1] * (transition[k - 1][0] *	// 0: match
		    f[i - 1][v] + transition[k - 1][4] * f[i - 1][v + 1] + transition[k - 1][7] * f[i - 1][v + 2]);
			
			set_u(w, bw, i - 1, k - ref_begin);
			if (w < bw2) f[i][u + 1] = emission[k][temp] * (transition[k][1] * f[i - 1][w] + transition[k][5] * f[i - 1][w + 1]);	// 1: insertion
	
			set_u(v, bw, i, k - 1 - ref_begin);
			if (i < read_len - 1 && k < window_len && v >= 0)
				f[i][u + 2] = transition[k - 1][2] * f[i][v] + transition[k - 1][8] * f[i][v + 2];	// 2: deletion
			
			s[i] += f[i][u] + f[i][u + 1] + f[i][u + 2];
		}

		/* rescale */
		for (k = beg; k <= end; k ++) {
			set_u(u, bw, i, k - ref_begin);
//fprintf(stderr, "s[%d]: %g\n", i, s[i]);
			f[i][u] /= s[i];	// 0: match
			f[i][u + 1] /= s[i];	// 1: insertion
			f[i][u + 2] /= s[i];	// 2: deletion
		}
	}

	/* sum of all forward path */
	f_final = 0;

	for (k = beg; k <= end; ++k) {
		set_u(u, bw, read_len - 1, k - ref_begin);
		if (u < 0 || u >= bw2) continue;
		f_final += transition[k][3] * f[read_len - 1][u] + transition[k][6] * f[read_len - 1][u + 1];
//fprintf(stderr, "t[%d][3]: %g\tf[%d][%d]: %g\tt[%d][6]: %g\tf[%d][%d]: %g\n", k, transition[k][3], read_len - 1, u, f[read_len - 1][u], k, transition[k][6], read_len - 1, u + 1, f[read_len - 1][u + 1]);
	}

	s[read_len] = f_final;

#ifdef VERBOSE_DEBUG
	f_final /= s[read_len];	// Debug
	fprintf(stderr, "f_final: %g\n", f_final);	// Debug
#endif

	/*--------------------*
	 * backword algorithm *
	 *--------------------*/
	x = ref_begin + read_len - 1 - bw; beg = x < 0 ? 0 : x; 
	x = ref_begin + read_len - 1 + bw; end = x < window_len ? x : window_len;
	/* b_0_E needs to be rescaled by s[read_len] */
	for (k = end; k >= beg; k --) {
		set_u(u, bw, read_len - 1, k - ref_begin);
		if (k > 0) b[read_len - 1][u] = transition[k][3] / s[read_len];	// 0: match
		b[read_len - 1][u + 1] = transition[k][6] / s[read_len];	// 1: insertion
	}

	/* rescale */
	for (k = end; k >= beg; k --) {
		set_u(u, bw, read_len - 1, k - ref_begin);
		b[read_len - 1][u] /= s[read_len - 1];	// 0: match
		b[read_len - 1][u + 1] /= s[read_len - 1];	// 1: insertion

#ifdef VERBOSE_DEBUG
		/* Debug: posterior probability */
		pp += b[read_len - 1][u] * f[read_len - 1][u] + b[read_len - 1][u + 1] * f[read_len - 1][u + 1];	// Debug
#endif

	}

#ifdef VERBOSE_DEBUG
	pp *= s[read_len - 1];	// Debug
	fprintf (stderr, "pp: %f\n", pp);	// Debug
#endif
	
	if (read_len > 1) {
		for (i = read_len - 2; i >= 0; i --) {
			int32_t even;
			temp1 = bam1_seqi(read, i + 1);
			temp = temp1 + pow(-1, temp1%2);
			x = ref_begin + i - bw; beg = 0 > x ? 0 : x;
			x = ref_begin + i + bw; end = x < window_len ? x : window_len; //	band end

			set_u(u, bw, i, end - ref_begin);
			set_u(v, bw, i + 1, end + 1 - ref_begin);
			set_u(w, bw, i + 1, end - ref_begin);
			set_u(y, bw, i, end + 1 - ref_begin);
			b[i][u] = transition[end][1] * emission[end][temp] * b[i + 1][w + 1];
			b[i][u + 1] = transition[end][5] * emission[end][temp] * b[i + 1][w + 1];
			if (end < window_len) {
				b[i][u] += transition[end][0] * emission[end + 1][temp1] * b[i + 1][v];
				b[i][u + 1] += transition[end][4] * emission[end + 1][temp1] * b[i + 1][v];	// 1: insertion
				b[i][u + 2] = transition[end][7] * emission[end + 1][temp1] * b[i + 1][v];
			}
			if (end < window_len - 1 && y < bw2) {
				b[i][u] += transition[end][2] * b[i][y + 2];	// 0: match
				b[i][u + 2] += transition[end][8] * b[i][y + 2];	// 2: deletion
			}

			even = beg > 2 ? beg : 2;
			for (k = end - 1; k >= even; k --) {
				set_u(u, bw, i, k - ref_begin);
				set_u(v, bw, i + 1, k + 1 - ref_begin);
				set_u(w, bw, i + 1, k - ref_begin);
				set_u(y, bw, i, k + 1 - ref_begin);
				if (w >= 0) {
					b[i][u] = transition[k][0] * emission[k + 1][temp1] * b[i + 1][v] +
					transition[k][1] * emission[k][temp] * b[i + 1][w + 1] + transition[k][2] * b[i][y + 2];	// 0: match

					b[i][u + 1] = transition[k][4] * emission[k + 1][temp1] *  
					b[i + 1][v] + transition[k][5] * emission[k][temp] * b[i + 1][w + 1];	// 1: insertion
				}
				else {
					b[i][u] = transition[k][0] * emission[k + 1][temp1] * b[i + 1][v] + transition[k][2] * b[i][y + 2];

					b[i][u + 1] = transition[k][4] * emission[k + 1][temp1] * b[i + 1][v];	// 1: insertion
				}

				b[i][u + 2] = transition[k][7] * emission[k + 1][temp1] *  
				b[i + 1][v] + transition[k][8] * b[i][y + 2];	// 2: deletion
			}

			if (beg <= 1) {
				set_u(u, bw, i, 1 - ref_begin);
				set_u(v, bw, i + 1, 2 - ref_begin);
				set_u(w, bw, i + 1, 1 - ref_begin);
				set_u(y, bw, i, 2 - ref_begin);
				if (w >= 0) {
					b[i][u] = transition[1][0] * emission[2][temp1] * b[i + 1][v] +
					transition[1][1] * emission[1][temp] * b[i + 1][w + 1] + transition[1][2] * b[i][y + 2];	// 0: match

					b[i][u + 1] = transition[1][4] * emission[2][temp1] *  
					b[i + 1][v] + transition[1][5] * emission[1][temp] * b[i + 1][w + 1];	// 1: insertion
				} else {
					b[i][u] = transition[1][0] * emission[2][temp1] * b[i + 1][v] + transition[1][2] * b[i][y + 2];	// 0: match

					b[i][u + 1] = transition[1][4] * emission[2][temp1] * b[i + 1][v];	// 1: insertion
				}
			}	

			if (beg == 0) {
				set_u(u, bw, i, 0 - ref_begin);
				set_u(v, bw, i + 1, 1 - ref_begin);
				set_u(w, bw, i + 1, 0 - ref_begin);
				if (w >= 0) b[i][u + 1] = transition[0][4] * emission[1][temp1] * b[i + 1][v] +
					transition[0][5] * emission[0][temp] * b[i + 1][w + 1];	// 1: insertion
				else b[i][u + 1] = transition[0][4] * emission[1][temp1] * b[i + 1][v];
			}

			/* rescale */
#ifdef VERBOSE_DEBUG
			pp = 0;	// Debug
#endif
 
			for (k = beg; k <= end; k ++) {
				set_u(u, bw, i, k - ref_begin);
				b[i][u] /= s[i];	// 0: match
				b[i][u + 1] /= s[i];	// 1: insertion
				b[i][u + 2] /= s[i];	// 2: deletion

#ifdef VERBOSE_DEBUG
				pp += b[i][u] * f[i][u] + b[i][u + 1] * f[i][u + 1]; // Debug
#endif
			}

#ifdef VERBOSE_DEBUG
			pp *= s[i];	// Debug
			fprintf (stderr, "pp: %f\ts[%d]: %g\n", pp, i, s[i]);	// Debug
#endif
		}

	}

	/* rescale */
	b_final = 0;
	temp1 = bam1_seqi(read, 0);
	temp = temp1 + pow(-1, temp1%2);

	for (k = beg; k <= end; ++k) {
		set_u(u, bw, 0, k - ref_begin);
		if (u < 0 || u >= bw2) continue;
		if (k > 0) b_final += transition[k - 1][9] * emission[k][temp1] * b[0][u] + 
		transition[k][10] * emission[k][temp] * b[0][u + 1];
		else b_final += transition[k][10] * emission[k][temp] * b[0][u + 1];
	}

#ifdef VERBOSE_DEBUG
	fprintf (stderr, "b_final: %g\n", b_final);	// Debug: b_final should equal to 1 

	// Debug: posterior probability for transition (each edge, no band) 
	{
		double pp_t = 0;
		for (i = 0; i < read_len - 1; i ++) {
			pp_t = 0;
			temp1 = bam1_seqi(read, i + 1);
			temp = temp1 + pow(-1, temp1%2);

			x = ref_begin + i - bw; beg = x > 0 ? x : 0;
			x = ref_begin + i + bw; end = window_len < x ? window_len : x; //	band end
			for (k = 2; k < window_len; k ++) {
				set_u(u, bw, i, k - ref_begin);
				set_u(v, bw, i + 1, k + 1 - ref_begin);
				if (u >= 0 && u < bw2 && v >= 0 && v < bw2)
					pp_t += f[i][u + 2] * transition[k][7] * emission[k + 1][bam1_seqi(read, i + 1)] * b[i + 1][v];
			}
			for (k = 1; k < window_len; k ++) {
				set_u(u, bw, i, k - ref_begin);
				set_u(v, bw, i + 1, k + 1 - ref_begin);
				if (u >= 0 && u < bw2 && v >= 0 && v < bw2) {
					pp_t += f[i][u] * transition[k][0] * emission[k + 1][bam1_seqi(read, i + 1)] * b[i + 1][v];
					pp_t += f[i][u + 1] * transition[k][4] * emission[k + 1][bam1_seqi(read, i + 1)] * b[i + 1][v];
				}

				set_u(v, bw, i + 1, k - ref_begin);
				if (u >= 0 && u < bw2 && v >= 0 && v < bw2) {
					pp_t += emission[k][temp] * f[i][u] * transition[k][1] * b[i + 1][v + 1];
					pp_t += emission[k][temp] * f[i][u + 1] * transition[k][5] * b[i + 1][v + 1];
				}
			}

			set_u(u, bw, i, 0 - ref_begin);
			set_u(v, bw, i + 1, 1 - ref_begin);
			if (u >= 0 && u < bw2 && v >= 0 && v < bw2)
				pp_t += f[i][u + 1] * transition[0][4] * emission[1][bam1_seqi(read, i + 1)] * b[i + 1][v];

			set_u(w, bw, i, window_len - ref_begin);
			set_u(y, bw, i + 1, window_len - ref_begin);
			if (w >= 0 && w < bw2 && y >= 0 && y < bw2) {
				pp_t += emission[window_len][temp] * f[i][w] * transition[window_len][1] * b[i + 1][y + 1];
				pp_t += emission[window_len][temp] * f[i][w + 1] * transition[window_len][5] * b[i + 1][y + 1];
			}			

			set_u(v, bw, i + 1, 0 - ref_begin);
			if (u >= 0 && u < bw2 && v >= 0 && v < bw2)
				pp_t += emission[0][temp] * f[i][u + 1] * transition[0][5] * b[i + 1][v + 1];
			fprintf (stderr, "i: %d\tpp_t: %g\n", i, pp_t);
		}
	}  // Debug
#endif
		
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

void baum_welch (double** transition, 
				 double** emission, 
				 int32_t window_begin,	// 0-based coordinate 
				 int32_t window_len, 
				 int32_t bw, 
				 reads* r, 
				 double df) {	// conversion condition: difference between 2 runs of EM

	double Pr = 10e100, diff = 1, p = 0;
	int32_t i, k, j, count = 0;
	double** t = calloc (window_len + 1, sizeof(double*));
	double** e = calloc (window_len + 1, sizeof(double*));
	for (i = 0; i <= window_len; i ++) {
		t[i] = calloc (16, sizeof(double));
		e[i] = calloc (16, sizeof(double));
	}		
	double** s_t = malloc ((window_len + 1) * sizeof(double*));
	for (k = 0; k < window_len + 1; k ++) {
		s_t[k] = malloc (3 * sizeof(double));
	}
	double** s_e = malloc ((window_len + 1) * sizeof(double*));
	for (k = 0; k < window_len + 1; k ++) {
		s_e[k] = malloc (2 * sizeof(double));
	}

	while (diff > df && count < 10) {
		if (count > 0 && p > 0) {
			for (k = 0; k <= window_len; k ++) {
				transition[k][0] = t[k][0];
				transition[k][1] = t[k][1]; 
				transition[k][2] = t[k][2];
				transition[k][4] = t[k][4];
				transition[k][5] = t[k][5];
				transition[k][7] = t[k][7];
				transition[k][8] = t[k][8];
				for (i = 0; i < 16; i ++) emission[k][i] = e[k][i];
			}
		}
		
		if (p > 0) {
			// Initialize new transition and emission matrixes. 
			for (k = 0; k <= window_len; k ++) {
				t[k][3] = transition[k][3];
				t[k][6] = transition[k][6];
				t[k][9] = transition[k][9];
				t[k][10] = transition[k][10];
				t[k][0] = t[k][1] = t[k][2] = t[k][4] = t[k][5] = t[k][7] = t[k][8] = 0; 
				for (i = 0; i < 16; i ++) {
					e[k][i] = 0; 
				}
			}
		}
		p = 0;	// likelihood 
		int32_t total_hl = 0;

		// Transition and emission matrixes training by a block of reads.
		for (j = 0; j < r->count; j ++) {
			uint8_t mq = r->qual[j];
fprintf(stderr, "count: %d\tmq: %d\n", count, mq);
fprintf(stderr, "j: %d\n", j);
			uint8_t* read_seq = &r->seqs[total_hl];
			total_hl += r->seq_l[j]/2 + r->seq_l[j]%2;
			if (count%(29/mq + 1) > 0) continue;
			int32_t read_len = r->seq_l[j];
			int32_t ref_begin = r->pos[j] + 1 - window_begin;
			int32_t bw2 = 3*(bw * 2 + 1);
			int32_t temp1, beg_i, end_i;
			double temp;
			double** f = (double**)calloc(read_len, sizeof(double*));
			double** b = (double**)calloc(read_len, sizeof(double*));
			for (i = 0; i < read_len; ++i) {
				f[i] = (double*)calloc(bw2, sizeof(double));
				b[i] = (double*)calloc(bw2, sizeof(double));
			}
			double* s = (double*)calloc(read_len + 1, sizeof(double));

			p += forward_backward (transition, emission, ref_begin, window_len, read_seq, read_len, f, b, s, bw);

			for (k = 0; k <= window_len; k ++) {
				beg_i = k - ref_begin - bw > 0 ? k - ref_begin - bw : 0;
				end_i = k - ref_begin + bw > read_len - 1 ? read_len - 1 : k - ref_begin + bw;
				for (i = beg_i; i < end_i; i ++) {
					int32_t u, v11, v01, v10;
					set_u(u, bw, i, k - ref_begin);
					set_u(v11, bw, i + 1, k + 1 - ref_begin);
					set_u(v01, bw, i, k + 1 - ref_begin);
					set_u(v10, bw, i + 1, k - ref_begin);

					temp1 = bam1_seqi(read_seq, i + 1);
					temp = emission[k][temp1 + (int32_t)pow(-1, temp1%2)];

					if (k < window_len) {
						t[k][0] += f[i][u] * transition[k][0] * emission[k + 1][temp1] * b[i + 1][v11];	// M_k -> M_k+1 
					
						t[k][4] += f[i][u + 1] * transition[k][4] * emission[k + 1][temp1] * b[i + 1][v11];	// I_k -> M_k+1 
					
						t[k][7] += f[i][u + 2] * transition[k][7] * emission[k + 1][temp1] * b[i + 1][v11];	// D_k -> M_k+1 
					
						if (v01 < bw2) {	
							t[k][2] += f[i][u] * transition[k][2] * b[i][v01 + 2] * s[i];	// M_k -> D_k+1 
			
							t[k][8] += f[i][u + 2] * transition[k][8] * b[i][v01 + 2] * s[i];	// D_k -> D_k+1
						}
					}					

					t[k][1] += f[i][u] * transition[k][1] * temp * b[i + 1][v10 + 1]; // M_k -> I_k 
			
					t[k][5] += f[i][u + 1] * transition[k][5] * temp * b[i + 1][v10 + 1];	// I_k -> I_k 
				}

				// i = read_len - 1 
				if (end_i > 1) {
					int32_t u, v01; 
					set_u(u, bw, end_i, k - ref_begin);
					set_u(v01, bw, end_i, k + 1 - ref_begin);
					if (v01 < bw2) {
						t[k][2] += f[end_i][u] * transition[k][2] * b[end_i][v01 + 2] * s[end_i];	// M_k -> D_k+1 
						t[k][8] += f[end_i][u + 2] * transition[k][8] * b[end_i][v01 + 2] * s[end_i];	// D_k -> D_k+1
					} 
				}
				
				for (i = beg_i; i <= end_i; i ++) {
					int32_t u, v01; 
					set_u(u, bw, i, k - ref_begin);
					set_u(v01, bw, i, k + 1 - ref_begin);

					temp1 = bam1_seqi(read_seq, i);
					if (k < window_len && v01 < bw2) e[k + 1][temp1] += f[i][v01] * b[i][v01] * s[i];	// M_k+1
					
					e[k][temp1 + (int32_t)pow(-1, temp1%2)] += f[i][u + 1] * b[i][u + 1] * s[i];	// I_k 
				}
			}

			beg_i = window_len - ref_begin - bw > 0 ? window_len - ref_begin - bw : 0;
			end_i = window_len - ref_begin + bw - 1 > read_len - 1 ? read_len - 1 : window_len - ref_begin + bw - 1;
			for (i = beg_i; i < end_i; i ++) {
				int32_t u, v10;
				set_u(u, bw, i, window_len - ref_begin);
				set_u(v10, bw, i + 1, window_len - ref_begin);

				temp1 = bam1_seqi(read_seq, i + 1);
				temp = emission[window_len][temp1 + (int32_t)pow(-1, temp1%2)];
				 // M_k -> I_k 
				t[window_len][1] += f[i][u] * transition[window_len][1] * temp * b[i + 1][v10 + 1];
					
				// I_k -> I_k 
				t[window_len][5] += f[i][u + 1] * transition[window_len][5] * temp * b[i + 1][v10 + 1];
			}
			for (i = beg_i; i <= end_i; i ++) {
				int32_t u;
				set_u(u, bw, i, window_len - ref_begin);
				temp1 = bam1_seqi(read_seq, i);
				e[window_len][temp1 + (int32_t)pow(-1, temp1%2)] += f[i][u + 1] * b[i][u + 1] * s[i];	// I_k 
			}
			free(s);
			for (i = 0; i < read_len; ++i) {
				free(f[i]);
				free(b[i]);
			}	
			free(f);
			free(b);
		}
		// Loop ending: Transition and emission matrixes training by a block of reads. 

		if (p > 0) {
			// Estimate transition probabilities. 
			s_t[0][1] = t[0][4] + t[0][5] + t[0][6];
			if (s_t[0][1] > 0) {
				t[0][4] /= s_t[0][1];
				t[0][5] /= s_t[0][1];
				t[0][6] /= s_t[0][1];
			}

			for (k = 1; k < window_len; k ++) {
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

	//fprintf(stderr, "k: %d\tt0: %g\tt1: %g\tt2: %g\tt3: %g\tt7: %g\tt8: %g\n", k, t[k][0], t[k][1], t[k][2], t[k][3], t[k][7], t[k][8]);
			}

			s_t[window_len][0] = t[window_len][1] + t[window_len][3];
			if (s_t[window_len][0] > 0) {
				t[window_len][1] /= s_t[window_len][0];
				t[window_len][3] /= s_t[window_len][0];
			}
				
			s_t[window_len][1] = t[window_len][5] + t[window_len][6];
			if (s_t[window_len][1] > 0) {
				t[window_len][5] /= s_t[window_len][1];
				t[window_len][6] /= s_t[window_len][1];
			}

			// Estimate emission probabilities. 
			// Match states 
			s_e[0][1] = e[0][0] + e[0][3] + e[0][5] + e[0][9] + e[0][14];
			if (s_e[0][1] > 0) {
				e[0][0] /= s_e[0][1];
				e[0][3] /= s_e[0][1];
				e[0][5] /= s_e[0][1];
				e[0][9] /= s_e[0][1];
				e[0][14] /= s_e[0][1];
			}
			for (k = 1; k <= window_len; k ++) {
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

			diff = fabs(Pr - p);
			Pr = p;
		}
		count ++;
	}
	for (k = 0; k <= window_len; k ++) {
		for (i = 0; i < 16; i ++) {
			transition[k][i] = t[k][i];
			emission[k][i] = e[k][i];
		//	fprintf(stderr, "t[%d][%d]: %g\t", k, i, transition[k][i]);
//			if (transition[k][2] > 0) fprintf(stderr, "t[%d][2]: %g\tref: %c\n", k, t[k][2], ref_seq[k - 1]);
		}
//		fprintf(stderr, "\n");
	}
	for (i = 0; i <= window_len; i ++) {
		free(s_e[i]);
		free(s_t[i]);
	    free(t[i]);
		free(e[i]);
	}
	free(s_e);
	free(s_t);
	free(t);		
	free(e);
} 
