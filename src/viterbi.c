/*
 * viterbi.c: Generate the insert sequence, genotype MNP and generate the best sequence alignment.
 * Author: Mengyao Zhao
 * Create date: 2012-05-17
 * Contact: zhangmp@bc.edu
 * Last revise: 2013-03-06 
 */

#include "bam.h"
#include "hmm.h"

#define set_u(u, b, i, k) (u)=((k)-(i)+(b))*3;

char* viterbi (double** transition, 
			   double** emission, 
			   double** v, 
			   int32_t ref_begin,
			   int32_t window_len,	// window size 
			   uint8_t* read, 
			   int32_t read_len,
			   int32_t bw, 
			   int32_t bpoint) {	// 1_based inside window (relative) coordinate

	int32_t i;	 /* iter of read 0_based */ 
	int32_t k;	 /* iter of reference 1_based */
	int32_t temp1, temp, u, w, x, bw2 = 3*(2*bw + 1); 
	int32_t beg = ref_begin - bw > 0 ? ref_begin - bw : 0, end = window_len < ref_begin + bw ? window_len : ref_begin + bw;
	double v_final, path1, path2;

	// v[0]
	temp1 = bam1_seqi(read, 0);
	temp = temp1 + pow(-1, temp1%2);
	set_u(u, bw, 0, beg - ref_begin);
	v[0][u + 1] = emission[beg][temp] * transition[beg][10];	// 1: insertion
	
	for (k = beg + 1; k <= end; ++k) {
		set_u(u, bw, 0, k - ref_begin);
		v[0][u] = emission[k][temp1] * transition[k - 1][9];	// 0: match
	}

	// v[i]
	for (i = 1; i < read_len; ++i) {
		int32_t r = ref_begin + i - bw;
		double max, path3;

		beg = r > 0 ? r : 0;
	
		// k = 0
		temp1 = bam1_seqi(read, i);
		temp = temp1 + pow(-1, temp1%2);
		set_u(u, bw, i, beg - ref_begin);
		v[i][u + 1] = emission[beg][temp] * v[i - 1][u + 1] * transition[beg][5];	// v[i,I_0]

		// k = 1
		set_u(w, bw, i, beg + 1 - ref_begin);
		v[i][w] = emission[beg + 1][temp1] * v[i - 1][u + 1] * transition[beg][4];	// v[i, M_1]
		
		path1 = v[i - 1][w] * transition[beg + 1][1];
		path2 = v[i - 1][w + 1] * transition[beg + 1][5];
		max = path1 > path2 ? path1 : path2;
		v[i][w + 1] = emission[beg + 1][temp] * max;	// v[i, I_1]

		// k = 2
		set_u(u, bw, i, beg + 2 - ref_begin);
		set_u(x, bw, i - 1, beg + 1 - ref_begin);
		path1 = v[i - 1][x] * transition[beg + 1][0];
		path2 = v[i - 1][x + 1] * transition[beg + 1][4];
		max = path1 > path2 ? path1 : path2;
		v[i][u] = emission[beg + 2][temp1] * max;	// v[i, M_2]

		set_u(x, bw, i - 1, beg + 2 - ref_begin);
		path1 = v[i - 1][x] * transition[beg + 2][1];
		path2 = v[i - 1][x + 1] * transition[beg + 2][5];
		max = path1 > path2 ? path1 : path2;
		v[i][u + 1] = emission[beg + 2][temp] * max;	// v[i, I_2]

		v[i][u + 2] = v[i][w] * transition[beg + 1][2];	// v[i, D_2]

		r = ref_begin + i + bw; end = window_len < r ? window_len : r; //	band end
		// k = 3 ... L - 1
		for (k = beg + 3; k < end; k ++) {
			set_u(u, bw, i, k - ref_begin);
			set_u(x, bw, i - 1, k - 1 - ref_begin);
			path1 = v[i - 1][x] * transition[k - 1][0];
			path2 = v[i - 1][x + 1] * transition[k - 1][4];
			path3 = v[i - 1][x + 2] * transition[k - 1][7];
			max = path1 > path2 ? path1 : path2;
			max = path3 > max ? path3 : max;
			v[i][u] = emission[k][temp1] * max;	// v[i, M_k]

			set_u(x, bw, i - 1, k - ref_begin);
			path1 = v[i - 1][x] * transition[k][1];
			path2 = v[i - 1][x + 1] * transition[k][5];
			max = path1 > path2 ? path1 : path2;
			v[i][u + 1] = emission[k][temp] * max;	// v[i, I_k]

			set_u(x, bw, i, k - 1 - ref_begin);
			path1 = v[i][x] * transition[k - 1][2];
			path2 = v[i][x + 2] * transition[k - 1][8];
			v[i][u + 2] = path1 > path2 ? path1 : path2;	// v[i, D_k]
		}

		// k = L
		set_u(u, bw, i, end - ref_begin);
		set_u(x, bw, i - 1, end - 1 - ref_begin);
		path1 = v[i - 1][x] * transition[end - 1][0];
		path2 = v[i - 1][x + 1] * transition[end - 1][4];
		path3 = v[i - 1][x + 2] * transition[end - 1][7];
		max = path1 > path2 ? path1 : path2;
		max = path3 > max ? path3 : max;
		v[i][u] = emission[end][temp1] * max;	// v[i, M_L]

		set_u(x, bw, i - 1, end - ref_begin);
		path1 = v[i - 1][x] * transition[end][1];
		path2 = v[i - 1][x + 1] * transition[end][5];
		max = path1 > path2 ? path1 : path2;
		v[i][u + 1] = emission[end][temp] * max;	// v[i, I_k]
	}
	
	//termination
	v_final = 0;
	for (k = 0; k <= window_len; ++k) {
		set_u(u, bw, read_len - 1, k - ref_begin);
		if (u < 0 || u >= bw2) continue;
		path1 = v[read_len - 1][u] * transition[k][3];
		path2 = v[read_len - 1][u + 1] * transition[k][6];
		v_final = path1 > v_final ? path1 : v_final;
		v_final = path2 > v_final ? path2 : v_final;
	}
}



