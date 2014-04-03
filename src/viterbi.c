/*
 * viterbi.c: Generate the insert sequence, genotype MNP and generate the best sequence alignment.
 * Author: Mengyao Zhao
 * Create date: 2012-05-17
 * Contact: zhangmp@bc.edu
 * Last revise: 2014-03-13
 */

#include <string.h>
#include <math.h>
#include <float.h>
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

typedef struct {
	int32_t* p;
	int32_t l;	// length of path
}p_path;

char num2base (int8_t num) {
	char base;
	switch (num) {
		case 0:
		case 1:
			base = 'A';
			break;
		case 2:
		case 3:
			base = 'C';
			break;
		case 4:
		case 5:
			base = 'G';
			break;
		case 8:
		case 9:
			base = 'T';
			break;
		case 14:
		case 15:
			base = 'N';
			break;
		default:
			fprintf(stderr, "The base number is assigned wrongly.\n");
			exit(1);
	} 
	return base;
}

p_path viterbi (double** transition, 
			   double** emission, 
			   int32_t ref_begin,	// relative read mapping location in window 1_based
			   int32_t window_len,	// window size 
			   uint8_t* read, 
			   int32_t read_len,
			   int32_t bw) { 

	int32_t i;	 /* iter of read 0_based */ 
	int32_t k;	 /* iter of reference 1_based */
	int32_t temp1, temp, u, w, x, bw2 = 3*(2*bw + 1), l = read_len, s_final, even; 
	int32_t beg = ref_begin - bw > 0 ? ref_begin - bw : 0, end = window_len < ref_begin + bw ? window_len : ref_begin + bw;
	p_path path;
	path.p = (int32_t*)malloc(l * sizeof(int32_t));
	double v_final, path1, path2;//, s;
	double** v = (double**)malloc(read_len * sizeof(double*));
	int32_t** state = (int32_t**)calloc(read_len, sizeof(int32_t*));

	for (i = 0; i < read_len; ++i) {
		v[i] = (double*)malloc(bw2 * sizeof(double));
		state[i] = (int32_t*)calloc(bw2, sizeof(int32_t));
		for (k = 0; k < bw2; ++k) {
			v[i][k] = -DBL_MAX;
			state[i][k] = -1;
		}
	}

	// v[0]
	temp1 = bam1_seqi(read, 0);
	temp = temp1 + pow(-1, temp1%2);

	// k = 0
	if (beg == 0) {
		set_u(u, bw, 0, 0 - ref_begin);
		v[0][u + 1] = log(emission[0][temp] * transition[0][10]);	// 1: insertion
	}
	
	// k = 1 ... L
	even = beg > 1 ? beg : 1;
	for (k = even; k <= end; ++k) {
		set_u(u, bw, 0, k - ref_begin);
		set_u(w, bw, 0, k - 1 - ref_begin);
		v[0][u] = log(emission[k][temp1] * transition[k - 1][9]);	// 0: match
		v[0][u + 1] = log(emission[k][temp] * transition[k][10]);	// 1: insertion
		if (k >= 2 && k < window_len) {
			path1 = v[0][w] + log(transition[k - 1][2]);
			path2 = v[0][w + 2] + log(transition[k - 1][8]);
			v[0][u + 2] = path1 > path2 ? path1 : path2;
			state[0][u + 2] = path1 > path2 ? w : w + 2;
		}	// v[0, D_k]
	}

	// v[i]
	for (i = 1; i < read_len; ++i) {
		int32_t r;
		double max, path3;
	
		r = ref_begin + i - bw; beg = r > 0 ? r : 0;
		temp1 = bam1_seqi(read, i);
		temp = temp1 + pow(-1, temp1%2);
		if (beg == 0) {	// k = 0
			set_u(u, bw, i, 0 - ref_begin);
			set_u(w, bw, i - 1, 0 - ref_begin);
			v[i][u + 1] = log(emission[0][temp] * transition[0][5]) + v[i - 1][w + 1];	// v[i,I_0]
			state[i][u + 1] = w + 1;
		} 

		if (beg <= 1) {	// k = 1
			set_u(w, bw, i, 1 - ref_begin);
			set_u(x, bw, i - 1, 0 - ref_begin);
			v[i][w] = log(emission[1][temp1] * transition[0][4]) + v[i - 1][x + 1];	// v[i, M_1]
			state[i][w] = x + 1;
			
			set_u(x, bw, i - 1, 1 - ref_begin);
			path1 = v[i - 1][x] + log(transition[1][1]);
			path2 = v[i - 1][x + 1] + log(transition[1][5]);
			max = path1 > path2 ? path1 : path2;
			v[i][w + 1] = log(emission[1][temp]) + max;	// v[i, I_1]
			state[i][w + 1] = path1 > path2 ? x : x + 1;
		}

		if (beg <= 2) {	// k = 2
			set_u(u, bw, i, 2 - ref_begin);
			set_u(w, bw, i, 1 - ref_begin);
			set_u(x, bw, i - 1, 1 - ref_begin);
			path1 = v[i - 1][x] + log(transition[1][0]);
			path2 = v[i - 1][x + 1] + log(transition[1][4]);
			max = path1 > path2 ? path1 : path2;
			v[i][u] = log(emission[2][temp1]) + max;	// v[i, M_2]

			set_u(x, bw, i - 1, 2 - ref_begin);
			if (x < bw2) {
				path1 = v[i - 1][x] + log(transition[2][1]);
				path2 = v[i - 1][x + 1] + log(transition[2][5]);
				max = path1 > path2 ? path1 : path2;
				v[i][u + 1] = log(emission[2][temp]) + max;	// v[i, I_2]
			} else v[i][u + 1] = log(emission[2][temp]);	// v[i, I_2]

			if (w >= 0) {
				v[i][u + 2] = v[i][w] + log(transition[1][2]);	// v[i, D_2]
				state[i][u + 2] = w;
			} 
		}

		r = ref_begin + i + bw; end = window_len < r ? window_len : r; //	band end
		// k = 3 ... L
		even = beg > 3 ? beg : 3;
		for (k = even; k <= end; k ++) {
			set_u(u, bw, i, k - ref_begin);
			set_u(x, bw, i - 1, k - 1 - ref_begin);
			path1 = v[i - 1][x] + log(transition[k - 1][0]);
			path2 = v[i - 1][x + 1] + log(transition[k - 1][4]);
			path3 = v[i - 1][x + 2] + log(transition[k - 1][7]);
			max = path1 > path2 ? path1 : path2;
			max = path3 > max ? path3 : max;
			v[i][u] = log(emission[k][temp1]) + max;	// v[i, M_k]
			state[i][u] = max == path1 ? x : (max == path2 ? x + 1 : x + 2);

			set_u(x, bw, i - 1, k - ref_begin);
			if (x < bw2) {
				path1 = v[i - 1][x] + log(transition[k][1]);
				path2 = v[i - 1][x + 1] + log(transition[k][5]);
				max = path1 > path2 ? path1 : path2;
				v[i][u + 1] = log(emission[k][temp]) + max;	// v[i, I_k]
				state[i][u + 1] = path1 > path2 ? x : x + 1;
			} 

			if (i < read_len - 1 && k < window_len) {
				set_u(x, bw, i, k - 1 - ref_begin);
				if (x >= 0) {
					path1 = v[i][x] + log(transition[k - 1][2]);
					path2 = v[i][x + 2] + log(transition[k - 1][8]);
					v[i][u + 2] = path1 > path2 ? path1 : path2;	// v[i, D_k]
					state[i][u + 2] = path1 > path2 ? x : x + 2;
				}
			}
		}
	}
	
	//termination
	v_final = -DBL_MAX;	// The smallest negative float number.
	s_final = -1; 
	for (k = 0; k <= window_len; ++k) {
		set_u(u, bw, read_len - 1, k - ref_begin);
		if (u < 0 || u >= bw2) continue;
			path1 = v[read_len - 1][u] + log(transition[k][3]);
			path2 = v[read_len - 1][u + 1] + log(transition[k][6]);
			s_final = path1 > v_final ? u : s_final; 
			v_final = path1 > v_final ? path1 : v_final;
			s_final = path2 > v_final ? u + 1 : s_final;
			v_final = path2 > v_final ? path2 : v_final;
	}

	// trace back
	u = s_final;
	x = 0;	// path index
	i = read_len - 1;
	while (k > 1 && i >= 0 && u > 0) {
		temp = u%3;
		temp1 = u/3;
		set_k(temp1, bw, i, k, ref_begin);
		if (x + 1 > l) {
			++l;
			kroundup32(l);
			path.p = realloc(path.p, l*sizeof(int32_t));
		}
		path.p[x++] = 3*k + temp;	// path is reversed
		u = state[i][u];
		if (temp == 0 || temp == 1) --i;
	}

	for (i = 0; i < read_len; ++i) {
		free(state[i]);
		free(v[i]);
	}
	free(state);
	free(v);

	path.l = x;
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
				kstring_t* ins,
				kstring_t* del,
				kstring_t* mva,
				khash_t(insert) *hi,	// key: 1-based relative position in window; value: insert_str1,insert_str2... (insert_str1 == insert_str2 is possible)
				khash_t(delet) *hd,
				khash_t(mnp) *hm) {

	int absent;
	khiter_t iter;
	
	if (ins && ins->l) {
		kputc(',', ins);
		iter = kh_put(insert, hi, pos, &absent);	// pos is a key, ret returns weather this key has existed
		if (absent == 1) {
			kh_value(hi, iter).l = kh_value(hi, iter).m = 0;
			kh_value(hi, iter).s = 0;
		}
		kputs(ins->s, &kh_value(hi, iter));
		ins->l = ins->m =0; ins->s = 0;
	} else if (mva && mva->l) {
		kputc(',', mva);
		iter = kh_put(mnp, hm, pos, &absent);	// pos is a key, ret returns weather this key has existed
		if (absent == 1) {
			kh_value(hm, iter).l = kh_value(hm, iter).m = 0;
			kh_value(hm, iter).s = 0;
		}
		kputs(mva->s, &kh_value(hm, iter));
		mva->l = mva->m = 0; mva->s = 0;
	} else if (del && del->l) {
		kputc(',', del);
		iter = kh_put(delet, hd, pos, &absent);	// pos is a key, ret returns weather this key has existed
		if (absent == 1) {
			kh_value(hd, iter).l = kh_value(hd, iter).m = 0;
			kh_value(hd, iter).s = 0;
		}
		kputs(del->s, &kh_value(hd, iter));
		del->l = del->m = 0; del->s = 0;
	}
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
		uint8_t* read_seq = &r->seqs[total_hl];
		total_hl += r->seq_l[j]/2 + r->seq_l[j]%2;
		int32_t c = 0;	// This is used to trace the read base.
		int32_t ref_begin = r->pos[j] + 1 - window_begin, i, k = 0, pos = 0, read_len = r->seq_l[j];
		p_path path = viterbi (transition, emission, ref_begin, window_len, read_seq, read_len, bw);

		kstring_t ins, del, mva;
		ins.l = ins.m = 0; ins.s = 0;
		del.l = del.m = 0; del.s = 0;
		mva.l = mva.m = 0; mva.s = 0;
		for (i = path.l - 1; i >= 0; --i) {	// Note: path is reversed
			int32_t read_base = bam1_seqi(read_seq, c);
			if (path.p[i] > 0 && path.p[i]%3)	{
				hash_seq (k, pos, 0, 0, &mva, hi, hd, hm);
				if (path.p[i]%3 == 1) {	// insert
					if (ins.l == 0) pos = path.p[i]/3;
					kputc(num2base(read_base), &ins);
					++c;
				} else {	// delet
					if (del.l == 0) pos = path.p[i]/3;
					kputc(ref_seq[path.p[i]/3 - 1], &del);
				}
			}else if (path.p[i]%3 == 0 && read_base != 15 && path.p[i]/3 > 0 && read_base != base2num(ref_seq, path.p[i]/3 - 1)) {	// mnp
				hash_seq (k, pos, &ins, &del, 0, hi, hd, hm);
				if (mva.l == 0) pos = path.p[i]/3;
				kputc(num2base(read_base), &mva);
				++c;
			} else if (path.p[i]%3 == 0 && path.p[i]/3 > 0 && read_base == base2num(ref_seq, path.p[i]/3 - 1)) {
				hash_seq (k, pos, &ins, &del, &mva, hi, hd, hm);	// Return to the main path
				++c;
			}
		}
		hash_seq (k, pos, &ins, &del, &mva, hi, hd, hm);
		free(path.p);

/*fprintf(stderr, "pos: %d\n", r->pos[j]);
	khiter_t d;	
	for (d = kh_begin(hd); d != kh_end(hd); ++d)
		if (kh_exist(hd, d)) fprintf(stderr, "key: %d\tvalue: %s\n", kh_key(hd, d), kh_value(hd, d).s);
	*/
	}

/*fprintf(stderr, "window_begin: %d\n", window_begin);
	khiter_t k;	
	for (k = kh_begin(hd); k != kh_end(hd); ++k)
		if (kh_exist(hd, k)) fprintf(stderr, "key: %d\tvalue: %s\n", kh_key(hd, k), kh_value(hd, k).s);
*/
}

