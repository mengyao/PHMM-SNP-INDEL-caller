/*
 * viterbi.c: Generate the insert sequence, genotype MNP and generate the best sequence alignment.
 * Author: Mengyao Zhao
 * Create date: 2012-05-17
 * Contact: zhangmp@bc.edu
 * Last revise: 2013-07-24
 */

#include <string.h>
#include <math.h>
#include "viterbi.h"
#include "bam.h"
#include "khash.h"
#include "kstring.h"

#define set_u(u, b, i, k) (u)=((k)-(i)+(b))*3;
#define set_k(u, b, i, k, r) {int x=(u)+(i)+(r)-(b); x=x>0?x:0; (k)=x;}
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

#ifndef KHASH
#define KHASH
KHASH_MAP_INIT_INT(insert, kstring_t)
KHASH_MAP_INIT_INT(mnp, kstring_t)
KHASH_MAP_INIT_INT(delet, kstring_t)
#endif

char num2base (int8_t num) {
	char base;
	switch (num) {
		case 1:
			base = 'A';
			break;
		case 2:
			base = 'C';
			break;
		case 4:
			base = 'G';
			break;
		case 8:
			base = 'T';
			break;
		case 15:
			base = 'N';
			break;
		default:
			fprintf(stderr, "The base number is assigned wrongly.\n");
			exit(1);
	} 
	return base;
}

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
	int32_t** state = (int32_t**)calloc(read_len, sizeof(int32_t*));

	for (i = 0; i < read_len; ++i) {
		v[i] = (double*)calloc(bw2, sizeof(double));
		state[i] = (int32_t*)calloc(bw2, sizeof(int32_t));
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
		set_u(w, bw, i - 1, beg - 1 - ref_begin);
		v[i][u + 1] = emission[beg][temp] * v[i - 1][w + 1] * transition[beg][5];	// v[i,I_0]
		state[i - 1][u + 1] = w + 1;

		// k = 1
		set_u(w, bw, i, beg + 1 - ref_begin);
		set_u(x, bw, i - 1, beg - ref_begin);
		v[i][w] = emission[beg + 1][temp1] * v[i - 1][x + 1] * transition[beg][4];	// v[i, M_1]
		state[i - 1][w] = x + 1;
		
		set_u(x, bw, i - 1, beg + 1 - ref_begin);
		path1 = v[i - 1][x] * transition[beg + 1][1];
		path2 = v[i - 1][x + 1] * transition[beg + 1][5];
		max = path1 > path2 ? path1 : path2;
		v[i][w + 1] = emission[beg + 1][temp] * max;	// v[i, I_1]
		state[i - 1][w + 1] = path1 > path2 ? x : x + 1;

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
//		fprintf(stderr, "i: %d, x: %d, beg: %d, bw: %d, bw2: %d, ref_begin: %d\n", i, x, beg, bw, bw2, ref_begin);
//		fprintf(stderr, "v[%d][%d]: %g\n", i - 1, x, v[i - 1][x]);
//		fprintf(stderr, "v[%d][%d]: %g\n", i - 1, x + 1, v[i - 1][x + 1]);
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
		if (x < bw2) {
			path1 = v[i - 1][x] * transition[end][1];
			max = path1 > path2 ? path1 : path2;
			v[i][u + 1] = emission[end][temp] * max;	// v[i, I_L]
			state[i - 1][u + 1] = path1 > path2 ? x : x + 1;
		} else {
			v[i][u + 1] = 0;
			state[i - 1][u + 1] = path1;
		}

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
		state[read_len - 1][0] = path1 > v_final ? u : state[read_len - 1][0]; 
		v_final = path1 > v_final ? path1 : v_final;
		state[read_len - 1][0] = path2 > v_final ? u + 1 : state[read_len - 1][0];
		v_final = path2 > v_final ? path2 : v_final;
	}

	// trace back
	temp = state[read_len - 1][0]%3;	
	temp1 = state[read_len - 1][0]/3;
	set_k(temp1, bw, read_len - 1, k, ref_begin);
	path[read_len - 1] = (3*k + temp);
	u = state[read_len - 1][0];
	for (i = read_len - 2; i >= 0; --i) {
		temp = state[i][u]%3;	// M: %3==0, I: %3==1, D: %3==2
		temp1 = state[i][u]/3;	// 1_based k= /3
		set_k(temp1, bw, i, k, ref_begin);
		fprintf(stderr, "pos: %d\tstate: %d\n", k, temp);
		path[i] = 3*k + temp;
		u = state[i][u];
	}

	for (i = 0; i < read_len; ++i) {
		free(state[i]);
		free(v[i]);
	}
	free(state);
	free(v);

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
			num = 0;
			break;
	}
	return num;
}

void hash_seq (int32_t k,
				int32_t pos,
				char* ins,
				char* del,
				char* mva,	
				khash_t(insert) *hi,	// key: 1-based relative position in window; value: insert_str1,insert_str2... (insert_str1 == insert_str2 is possible)
				khash_t(delet) *hd,
				khash_t(mnp) *hm) {

	int absent;
	khiter_t iter;
	
	if (ins) {
		ins[k++] = ',';
		ins[k] = '\0';
		iter = kh_put(insert, hi, pos, &absent);	// pos is a key, ret returns weather this key has existed
		if (absent == 1) {
			kh_value(hi, iter).l = kh_value(hi, iter).m = 0;
			kh_value(hi, iter).s = 0;
		}
		kputs(ins, &kh_value(hi, iter));
		free(ins);
		ins = 0;
	} else if (mva) {
		mva[k++] = ',';
		mva[k] = '\0';
		iter = kh_put(mnp, hm, pos, &absent);	// pos is a key, ret returns weather this key has existed
		if (absent == 1) {
			kh_value(hm, iter).l = kh_value(hm, iter).m = 0;
			kh_value(hm, iter).s = 0;
		}
		kputs(mva, &kh_value(hm, iter));
		free(mva);
		mva = 0;
	} else if (del) {
//	fprintf(stderr, "hash del pos: %d\tdel: %s\n", pos, del);
		del[k++] = ',';
		del[k] = '\0';
		iter = kh_put(delet, hd, pos, &absent);	// pos is a key, ret returns weather this key has existed
		if (absent == 1) {
			kh_value(hd, iter).l = kh_value(hd, iter).m = 0;
			kh_value(hd, iter).s = 0;
		}
		kputs(del, &kh_value(hd, iter));
		free(del);
		del = 0;
	}
}

// Update the indel or mnp string
char* update (char* str, char base, int32_t* max_len, int32_t* k, int32_t* pos, int32_t begin_pos) {
	if (! str) {
		*max_len = 128;
		str = malloc ((*max_len)*sizeof(char));
		*k = 0;
		*pos = begin_pos;	// 1_based k
	}
	if (*k + 2 > *max_len) { 
		++ (*max_len);
		kroundup32(*max_len);
		str = realloc(str, (*max_len) * sizeof(char));	
	}
	str[*k] = base;
	++ (*k);
//	fprintf(stderr, "str[*k - 1]: %c\n", str[*k - 1]);
//	fprintf(stderr, "str: %s\n", str);
	return str;
}

// Generate the hash of insertion, mnp and deletion
void hash_imd (double** transition, 
				double** emission, 				 
				char* ref_seq, 
				 int32_t window_begin,	// 0-based coordinate 
				 int32_t window_len, 
				 int32_t bw,
				 reads* r,
				khash_t(insert) *hi,	// key: 1-based relative position in window; value: insert_str1,insert_str2... (insert_str1 == insert_str2 is possible)
				khash_t(mnp) *hm,
				khash_t(delet) *hd) {

	int32_t j, total_hl = 0;
	for (j = 0; j < r->count; j ++) {
		char* ins = 0, *del = 0, *mva = 0;	// insertion, deletion and mnp seq
		uint8_t* read_seq = &r->seqs[total_hl];
		total_hl += r->seq_l[j]/2 + r->seq_l[j]%2;
		int32_t len_del = 128, len_mva = 128, len_ins = 128, i, k = 0, pos = 0, read_len = r->seq_l[j];
		int32_t ref_begin = r->pos[j] + 1 - window_begin;
		int32_t* path = viterbi (transition, emission, ref_begin, window_len, read_seq, read_len, bw);
		for (i = 0; i < read_len; ++i) {
			int32_t read_base = bam1_seqi(read_seq, i);
			if (path[i]%3)	{
				hash_seq (k, pos, 0, 0, mva, hi, hd, hm);
				if (path[i]%3 == 1) update(ins, num2base(read_base), &len_ins, &k, &pos, path[i]/3);	// insert
				else {
				//	fprintf(stderr, "ref[%d]: %c\n", path[i]/3 - 1, ref_seq[path[i]/3 - 1]);
					del = update(del, ref_seq[path[i]/3 - 1], &len_del, &k, &pos, path[i]/3);	// delet
				//	fprintf(stderr, "del: %s\n", del);
				}
			}else if (path[i]%3 == 0 && read_base != 15 && path[i]/3 > 0 && read_base != base2num(ref_seq, path[i]/3 - 1)) {	// mnp
				hash_seq (k, pos, ins, del, 0, hi, hd, hm);
				update(mva, num2base(read_base), &len_mva, &k, &pos, path[i]/3);
			} else if (path[i]%3 == 0 && path[i]/3 > 0 && read_base == base2num(ref_seq, path[i]/3 - 1)) {
			//	fprintf(stderr, "return, ins: %s, mva: %s, del: %s\n", ins, mva, del);
				hash_seq (k, pos, ins, del, mva, hi, hd, hm);	// Return to the main path
			}
		}
		hash_seq (k, pos, ins, del, mva, hi, hd, hm);
		free(path);
	}

	khiter_t k;	
	for (k = kh_begin(hd); k != kh_end(hd); ++k)
	//	if (kh_exist(hd, k)) fprintf(stderr, "%s\n", kh_value(hd, k).s);
}

