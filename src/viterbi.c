/*
 * viterbi.c: Generate the insert sequence, genotype MNP and generate the best sequence alignment.
 * Author: Mengyao Zhao
 * Create date: 2012-05-17
 * Contact: zhangmp@bc.edu
 * Last revise: 2013-03-21 
 */

#include <string.h>
#include "bam.h"
#include "hmm.h"
#include "khash.h"

#define set_u(u, b, i, k) (u)=((k)-(i)+(b)i)*3;
#define set_k(u, b, i, k) (k)=(u)/3+(i)-(b);

int32_t* viterbi (double** transition, 
			   double** emission, 
			   int32_t ref_begin,	// relative read mapping location in window 1_based
			   int32_t window_len,	// window size 
			   uint8_t* read, 
			   int32_t read_len,
			   int32_t bw) { 

	int32_t i;	 /* iter of read 0_based */ 
	int32_t k;	 /* iter of reference 1_based */
	int32_t temp1, temp, u, w, x, bw2 = 3*(2*bw + 1); 
	int32_t beg = ref_begin - bw > 0 ? ref_begin - bw : 0, end = window_len < ref_begin + bw ? window_len : ref_begin + bw;
	int32_t* path = (int32_t*)malloc(read_len * sizeof(int32_t));
	double v_final, path1, path2, s;
	double** v = (double**)calloc(read_len, sizeof(double*));
	double** state = (double**)calloc(read_len, sizeof(double*));

	for (i = 0; i < read_len; ++i) {
		v[i] = (double*)calloc(bw2, sizeof(double));
		state[i] = (double*)calloc(bw2, sizeof(double));
	}

	// v[0]
	temp1 = bam1_seqi(read, 0);
	temp = temp1 + pow(-1, temp1%2);

	// k = 0
	set_u(u, bw, 0, beg - ref_begin);
	v[0][u + 1] = emission[beg][temp] * transition[beg][10];	// 1: insertion
	s = v[0][u + 1];
	
	// k = 1 ... L
	for (k = beg + 1; k <= end; ++k) {
		set_u(u, bw, 0, k - ref_begin);
		v[0][u] = emission[k][temp1] * transition[k - 1][9];	// 0: match
		v[0][u + 1] = emission[k][temp] * transition[k][10];	// 1: insertion
		s += v[0][u] + v[0][u + 1];
	}

	/* rescale */
	for (k = beg; k <= end; k ++) {
		set_u(u, bw, 0, k - ref_begin);
		v[0][u] /= s;	// 0: match
		v[0][u + 1] /= s;	// 1: insertion
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
		state[i - 1][u + 1] = u + 1;


		// k = 1
		set_u(w, bw, i, beg + 1 - ref_begin);
		v[i][w] = emission[beg + 1][temp1] * v[i - 1][u + 1] * transition[beg][4];	// v[i, M_1]
		state[i - 1][w] = u + 1;
		
		path1 = v[i - 1][w] * transition[beg + 1][1];
		path2 = v[i - 1][w + 1] * transition[beg + 1][5];
		max = path1 > path2 ? path1 : path2;
		v[i][w + 1] = emission[beg + 1][temp] * max;	// v[i, I_1]
		state[i - 1][w + 1] = path1 > path2 ? w : w + 1;

		s = v[i][u + 1] + v[i][w] + v[i][w + 1];

		// k = 2
		set_u(u, bw, i, beg + 2 - ref_begin);
		set_u(x, bw, i - 1, beg + 1 - ref_begin);
		path1 = v[i - 1][x] * transition[beg + 1][0];
		path2 = v[i - 1][x + 1] * transition[beg + 1][4];
		max = path1 > path2 ? path1 : path2;
		v[i][u] = emission[beg + 2][temp1] * max;	// v[i, M_2]
		state[i - 1][u] = path1 > path2 ? x : x + 1;

		set_u(x, bw, i - 1, beg + 2 - ref_begin);
		path1 = v[i - 1][x] * transition[beg + 2][1];
		path2 = v[i - 1][x + 1] * transition[beg + 2][5];
		max = path1 > path2 ? path1 : path2;
		v[i][u + 1] = emission[beg + 2][temp] * max;	// v[i, I_2]
		state[i - 1][u + 1] = path1 > path2 ? x : x + 1;

		v[i][u + 2] = v[i][w] * transition[beg + 1][2];	// v[i, D_2]
		state[i - 1][u + 2] = w;

		s += v[i][u] + v[i][u + 1] + v[i][u + 2];

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
			state[i - 1][u] = max == path1 ? x : (max == path2 ? x + 1 : x + 2);

			set_u(x, bw, i - 1, k - ref_begin);
			path1 = v[i - 1][x] * transition[k][1];
			path2 = v[i - 1][x + 1] * transition[k][5];
			max = path1 > path2 ? path1 : path2;
			v[i][u + 1] = emission[k][temp] * max;	// v[i, I_k]
			state[i - 1][u + 1] = path1 > path2 ? x : x + 1;

			set_u(x, bw, i, k - 1 - ref_begin);
			path1 = v[i][x] * transition[k - 1][2];
			path2 = v[i][x + 2] * transition[k - 1][8];
			v[i][u + 2] = path1 > path2 ? path1 : path2;	// v[i, D_k]
			state[i - 1][u + 2] = path1 > path2 ? x : x + 2;
	
			s += v[i][u] + v[i][u + 1] + v[i][u + 2];
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
		state[i - 1][u] = max == path1 ? x : (max == path2 ? x + 1 : x + 2);

		set_u(x, bw, i - 1, end - ref_begin);
		path1 = v[i - 1][x] * transition[end][1];
		path2 = v[i - 1][x + 1] * transition[end][5];
		max = path1 > path2 ? path1 : path2;
		v[i][u + 1] = emission[end][temp] * max;	// v[i, I_k]
		state[i - 1][u + 1] = path1 > path2 ? x : x + 1;

		s += v[i][u] + v[i][u + 1];

		/* rescale */
		for (k = beg; k <= end; k ++) {
			set_u(u, bw, i, k - ref_begin);
			v[i][u] /= s;	// 0: match
			v[i][u + 1] /= s;	// 1: insertion
			v[i][u + 2] /= s;	// 2: deletion
		}
	}
	
	//termination
	v_final = 0;
	for (k = 0; k <= window_len; ++k) {
		set_u(u, bw, read_len - 1, k - ref_begin);
		if (u < 0 || u >= bw2) continue;
		path1 = v[read_len - 1][u] * transition[k][3];
		path2 = v[read_len - 1][u + 1] * transition[k][6];
		v_final = path1 > v_final ? path1 : v_final;
		state[read_len - 1][0] = path1 > v_final ? u : state[read - 1][0]; 
		v_final = path2 > v_final ? path2 : v_final;
		state[read_len - 1][0] = path2 > v_final ? u + 1 : state[read - 1][0];
	}

	for (i = 0; i < read_len; ++i) {
		free(state[i]);
		free(v[i]);
	}
	free(state);
	free(v);

	// trace back
	temp = state[read - 1][0]%3;
	temp1 = state[read - 1][0]/3;
	set_k(temp1, bw, read_len - 1, k);
	path[read_len - 1] = 3*k + temp;
	u = state[read - 1][0];
	for (i = read_len - 2; i >= 0; --i) {
		temp = state[i][u]%3;	// M: %3==0, I: %3==1, D: %3==2
		temp1 = state[i][u]/3;	// 1_based k= /3
		set_k(temp1, bw, i, k);
		path[i] = 3*k + temp;
		u = state[i][u];
	}

	return path;
}

int32_t base2num (char* seq, int32_t k) {
	int32_t num;
	switch (seq[k]) {
		case 'A':
		case 'a':
			num = 1;
			break;
		case 'C':
		case 'c':
			num = 2;
			break;
		case 'G':
		case 'g':
			num = 4;
			break;
		case 'T':
		case 't':
			num = 8;
			break;
		default:
			fprintf(stderr, "Wrong reference sequence. \n");
			exit (1);
			break;
	}
	return num;
}

void hash_insert_mnp (double** transition, 
				double** emission, 				 
				char* ref_seq, 
				 int32_t window_begin,	// 0-based coordinate 
				 int32_t window_len, 
				 int32_t bw, 
				 reads* r,
				khash_t(insert) *hi,
				khash_t(mnp) *hm) {

	int32_t j, total_hl = 0;
	int ret;
	khiter_t k;
	for (j = 0; j < r->count; j ++) {
		uint8_t* read_seq = &r->seqs[total_hl];
		total_hl += r->seq_l[j]/2 + r->seq_l[j]%2;
		int32_t i, k, pos, read_len = r->seq_l[j], flag = 0;	// flag == 0: no variation, flag == 1: insertion, flag == 2: mnp
		int32_t ref_begin = r->pos[j] + 1 - window_begin;
		int32_t* path = viterbi (transition, emission, ref_begin, window_len, read_seq, read_len, bw);
		for (i = 0; i < read_len; ++i) {
			int32_t read_base = bam1_seqi(read_seq, i);
			if (path[i]%3 == 1)	{	// insert
				if (flag == 0) {
					char* var = malloc ((read_len + 2) * sizeof(char));
					k = 0;
					var[k++] = ',';
					pos = path[i]/3;	// 1_based k
					flag = 1;
				} else if (flag == 2) {
					var[k] = '\0';
					k = kh_put(mnp, hm, pos, &ret);	// pos is a key, ret returns weather this key has existed
					if (ret == 1) {
						k = kh_get(mnp, hm, pos); 
						kh_value(hm, k) = strcat(kh_value(hm, k), var);
					} else kh_value(hm, k) = var;
					free(var);		
					char* var = malloc ((read_len + 2) * sizeof(char));
					k = 0;
					var[k++] = ',';
					pos = path[i]/3;	// 1_based k
					flag = 1;
				}
				var[k++] = read_base;
			}else if (path[i]%3 == 0 && read_base != base2num(ref_seq[path[i]/3 - 1])) {	// mnp
				if (flag == 0) {
					char* var = malloc ((read_len + 2) * sizeof(char));
					k = 0;
					var[k++] = ',';
					pos = path[i]/3;	// 1_based k
					flag = 2;
				} else if (flag == 1) {
					var[k] = '\0';
					k = kh_put(insert, hi, pos, &ret);	// pos is a key, ret returns weather this key has existed
					if (ret == 1) {
						k = kh_get(insert, hi, pos); 
						kh_value(hi, k) = strcat(kh_value(hi, k), var);
					} else kh_value(hi, k) = var;
					free(var);		
					char* var = malloc ((read_len + 2) * sizeof(char));
					k = 0;
					var[k++] = ',';
					pos = path[i]/3;	// 1_based k
					flag = 2;
				}
				var[k++] = read_base;
			}
		}
		free(path);
	}
}

