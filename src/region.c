/*
 * region.c: Get reference and alignments in a region using samtools-0.1.18
 * Author: Mengyao Zhao
 * Create date: 2011-06-05
 * Last revise date: 2014-08-11
 * Contact: zhangmp@bc.edu 
 */

#include <zlib.h>
#include <stdlib.h>
#include <string.h>
#include "khash.h"
#include "hmm.h"
#include "sicall.h"
#include "viterbi.h"
#include "region.h"
#include "pileup.h"

/*! @function
  @abstract  Round an integer to the next closest power-2 integer.
  @param  x  integer to be rounded (in place)
  @discussion x will be modified.
 */
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#define WINDOW_EDGE 50
//#define WINDOW_EDGE 20
#define WINDOW_SIZE 1000

/* This table is used to transform nucleotide letters into numbers. */
const int8_t nt_table[128] = {
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

typedef struct {
	double** transition;
	double** emission;
} profile;

int32_t buffer_read1 (bam1_t* bam, reads* r, int32_t window_begin, int32_t window_end, int32_t* count, int32_t* half_len) {
	int32_t read_len = bam->core.l_qseq, char_len, j;
	uint32_t* cigar = bam1_cigar(bam);
	uint8_t* read_seq = bam1_seq(bam);

	// Buffer reads in this region.
	if (bam->core.pos < window_begin) {
	// read head is aligned out of the window: truncate the head
		int32_t pos = bam->core.pos;
		int32_t clip_len = 0;
		uint16_t cigar_count = 0;
		while (pos < window_begin && cigar_count < bam->core.n_cigar) {
			int32_t operation = 0xf & *cigar;
			int32_t length;
			if (operation == 0 || operation == 7 || operation == 8) {	// M, =, X
				length = (0xfffffff0 & *cigar)>>4;
				if ((pos + length - 1) >= window_begin) {
					clip_len += window_begin - pos;
				} else clip_len += length;
				pos += length;
			} else if (operation == 1 || operation == 4) {	// I, S
				length = (0xfffffff0 & *cigar)>>4;
				clip_len += length;
			} else if (operation == 2 || operation == 3) {	// D, N
				length = (0xfffffff0 & *cigar)>>4;
				pos += length;
			}
			++ cigar;
			++ cigar_count;
		}
		if (clip_len%2) {	// Remove one more read sequence residual.
			read_len -= (clip_len + 1);
			if (read_len < 13) return 0;
			read_seq += ((clip_len + 1)/2);
			r->pos[*count] = window_begin + 1;
		} else {
			read_len -= clip_len;
			if (read_len < 13) return 0;
			read_seq += (clip_len/2);
			r->pos[*count] = window_begin;
		}
	} else if (bam->core.pos + read_len > window_end) {	// It will not meet this condition, when called by slide_window_whole	
		int32_t pos = bam->core.pos;
		uint16_t cigar_count = 0;
		read_len = 0;
		while (pos <= window_end && cigar_count < bam->core.n_cigar) {
			int32_t operation = 0xf & *cigar;
			int32_t length = 0;
			if (operation == 0 || operation == 7 || operation == 8) {	// M, =, X
				length = (0xfffffff0 & *cigar)>>4;
				if ((pos + length) > window_end) read_len += (window_end - pos + 1);
				else read_len += length;
				pos += length;
			} else if (operation == 1 || operation == 4) {	// I, S
				length = (0xfffffff0 & *cigar)>>4;
				read_len += length;
			} else if (operation == 2 || operation == 3) {	// D, N
				length = (0xfffffff0 & *cigar)>>4;
				pos += length;
			}
			++ cigar;
			++ cigar_count;
		}
		if (read_len < 13) return 0;
	}	
	if (bam->core.pos >= window_begin) r->pos[*count] = bam->core.pos;
	r->seq_l[*count] = read_len;
	char_len = read_len/2;
	for (j = *half_len; j < *half_len + char_len; j ++) {
		r->seqs[j] = read_seq[j - *half_len];
	}
	r->qual[*count] = bam->core.qual;
	if (read_len%2) r->seqs[j] = read_seq[j - *half_len];
	(*half_len) += char_len;
	if (read_len%2) (*half_len) ++;
	(*count) ++;
	
	return 1;
}

// Change a string to number
int32_t value(int8_t* num_seq, int32_t p, int32_t seg_len) {
	int32_t i, s = num_seq[p];
	for (i = p + 1; i < p + seg_len; ++i) s = (s << 2) | num_seq[i];
	return s;
}

// Mark the repeat region, and the repeated segment beginning position.
void insert_group (char* ref_seq, int32_t ref_len, profile* hmm, int32_t beg) {
	int32_t length = ref_len - beg + 1;
	int8_t* r_mark = (int8_t*)calloc (length, sizeof(int8_t));
	int8_t* num_seq = (int8_t*)malloc (sizeof(int8_t)*length);
	int32_t* v_s = (int32_t*)calloc (length, sizeof(int32_t)); 	// value of the short strings
	int32_t i, seg_len, pos_i = 0;
	double sum, max_i = 0.05, sum_i = 0;

	for (i = 0; i < length; ++i) {
		num_seq[i] = nt_table[(int)ref_seq[i + beg - 1]];
	}
	for (seg_len = 1; seg_len < 7; ++ seg_len) {
		memset (v_s, -1, length * sizeof(int32_t));
		for  (i = seg_len - 1; i < length - seg_len; ++i) {
			v_s[i] = value(num_seq, i, seg_len);
			if (i - 2*seg_len + 1 >= 0 && v_s[i] == v_s[i - seg_len]) r_mark[i] = r_mark[i - seg_len] = seg_len;
		}
	}

	free (v_s);
	free (num_seq);
	
	// Group insertion signal.
	for (i = 0; i < length; ++i) {
		int32_t jump = r_mark[i];
		// Signal group when tandem repeat / homopolymer region end.
		if (i > 1 && r_mark[i - 1] > 0 && jump != r_mark[i - 1]) {
			hmm->transition[pos_i][1] += sum_i;
			hmm->transition[pos_i][0] -= sum_i;
			hmm->transition[pos_i][0] = hmm->transition[pos_i][0] > 0 ? hmm->transition[pos_i][0] : 0;
			sum = hmm->transition[pos_i][0] + hmm->transition[pos_i][1] + hmm->transition[pos_i][2] + hmm->transition[pos_i][3];
			hmm->transition[pos_i][0] /= sum;
			hmm->transition[pos_i][1] /= sum;
			hmm->transition[pos_i][2] /= sum;
			hmm->transition[pos_i][3] /= sum;
			sum_i = 0;
		}

		// In tandem repeat / homopolymer region: gether the signal and seek the strongest signal.
		if (jump > 0 && hmm->transition[i + beg][1] > 0.05) {
			if (hmm->transition[i + beg][1] > max_i) {
				pos_i = i + beg;
				max_i = hmm->transition[i + beg][1];
			} 
			sum_i += hmm->transition[i + beg][1];
			hmm->transition[i + beg][0] += hmm->transition[i + beg][1];
			hmm->transition[i + beg][1] = 0;
		
			sum = hmm->transition[i + beg][0] + hmm->transition[i + beg][1] + hmm->transition[i + beg][2] + hmm->transition[i + beg][3];	
			hmm->transition[i + beg][0] /= sum;
			hmm->transition[i + beg][1] /= sum;
			hmm->transition[i + beg][2] /= sum;
			hmm->transition[i + beg][3] /= sum;
		}
	}	

	free (r_mark);
}

void call_var (bam_header_t* header,
				reads* r, 
				p_info* cinfo,
				char* ref_seq,
				int32_t tid,	
			   	  int32_t window_begin,	// 0-based 
			   	  int32_t window_end,
			   	  int32_t region_begin,	// -1: slide_window_whole
			   	  int32_t region_end,	// only used in slide_window_region 
			   	  int32_t size) {

	int32_t frame_begin, frame_end, temp, i, region_len, ref_len = strlen(ref_seq);
//	char* ref_seq = faidx_fetch_seq(fai, header->target_name[tid], window_begin, window_end, &ref_len);
	double** e = (double**)calloc(ref_len + size + 1, sizeof(double*));
	profile* hmm = (profile*)malloc(sizeof(profile));
	khash_t(insert) *hi = kh_init(insert);
	khash_t(mnp) *hm = kh_init(mnp);
	khash_t(delet) *hd = kh_init(delet);
	khiter_t k;

	if (region_end == 2147483647 || region_end == 536870912) region_end = window_begin + ref_len;	// slid_window_whole || slid_window_region user only gave the chromosome number
	region_len = region_end - region_begin;

	if (size > 2) hmm->transition = transition_init (0.1, 0.3, 0.2, 2.5, 0.4, ref_len);
	else hmm->transition = transition_init (0.01, 0.3, 0.2, 2.5, 0.4, ref_len);
	hmm->emission = emission_init(ref_seq, 0.24, 0.9, 0.32); 

	//Copy the initiated emission matrix for the Viterbi.
	for (k = 0; k <= ref_len; ++k) { 
		e[k] = (double*)calloc(16, sizeof(double));
		for (i = 0; i < 16; ++i)
			e[k][i] = hmm->emission[k][i];
	}

	baum_welch (hmm->transition, hmm->emission, window_begin, ref_len, size, r, 0.01);

	// Group the homopolymer deletion signal to the most left position.
	for (i = 0; i < ref_len - 3; ++i) {

		// Group homopolymer deletion signal.
		if (ref_seq[i] == ref_seq[i + 1] && ref_seq[i] == ref_seq[i + 2]) {
			double sum;
			int32_t j = i + 1;//, pos_i = j;
			while (ref_seq[j] == ref_seq[i]) {
				if (hmm->transition[j][2] > 0.001) { 	// Group deletion signal.
					hmm->transition[i][2] += hmm->transition[j][2];
					hmm->transition[j][0] += hmm->transition[j][2] - 0.001;
					hmm->transition[j][2] = 0.001;
					sum = hmm->transition[j][0] + hmm->transition[j][1] + hmm->transition[j][2] + hmm->transition[j][3];	
					hmm->transition[j][0] /= sum;
					hmm->transition[j][1] /= sum;
					hmm->transition[j][2] /= sum;
					hmm->transition[j][3] /= sum;
				}
				++j;
			}

			// Group deletion signal.
			sum = hmm->transition[i][0] + hmm->transition[i][1] + hmm->transition[i][2] + hmm->transition[i][3];
			hmm->transition[i][0] /= sum;
			hmm->transition[i][1] /= sum;
			hmm->transition[i][2] /= sum;
			hmm->transition[i][3] /= sum;
		}
	}
	
	i = 1;
	while (i < ref_len - 3 && hmm->transition[i][1] < 0.05) ++i;
	if (i < ref_len - 3) insert_group (ref_seq, ref_len, hmm, i - 1);
/*
	for (k = 0; k <= ref_len; ++k) {
		for (i = 0; i < 10; ++i) fprintf(stderr, "t[%d][%d]: %g\t", k, i, hmm->transition[k][i]);
		fprintf(stderr, "\n");
	}
*/
	hash_imd (hmm->transition, e, ref_seq, window_begin, ref_len, size, r, hi, hm, hd);

	if (region_begin >= 0 && region_len < WINDOW_SIZE) {	// small region
		frame_begin = region_begin;
		frame_end = region_end;
	} else { 
		temp = window_begin + WINDOW_EDGE;
		frame_begin = temp > region_begin ? temp : region_begin;
		temp = window_begin + ref_len - WINDOW_EDGE;
		frame_end = temp < region_end ? temp : region_end;
	}

	if(frame_end > frame_begin) {
		likelihood (header, hmm->transition, hmm->emission, ref_seq, cinfo, tid, window_begin, frame_begin, frame_end, size, 0, hi, hm, hd);
	}	

	for (k = kh_begin(hm); k != kh_end(hm); ++k)
		if (kh_exist(hm, k)) free(kh_value(hm, k).s);
	kh_destroy(mnp, hm);    		
	for (k = kh_begin(hi); k != kh_end(hi); ++k)
		if (kh_exist(hi, k)) free(kh_value(hi, k).s);
	kh_destroy(insert, hi);
	for (k = kh_begin(hd); k != kh_end(hd); ++k)
		if (kh_exist(hd, k)) free(kh_value(hd, k).s);
	kh_destroy(delet, hd);
    		
	transition_destroy(hmm->transition, ref_len);
	emission_destroy(hmm->emission, ref_len);
	free(hmm);
//	free(ref_seq);

	return;
}

// Use the new read to update the depth array and the summary of quality array. 
p_info* add_depth (p_info* cinfo, int32_t* d, int32_t read_beg, int32_t read_length, uint8_t qual) {
	int32_t i, beg;
	if(read_beg + read_length > *d) {
		int32_t orig = *d;
		(*d) = read_beg + read_length + 1;
		kroundup32(*d);
		cinfo = realloc(cinfo, (*d)*sizeof(p_info));
		for (i = orig; i < (*d); ++ i) {
			cinfo[i].depth = 0;
			cinfo[i].mqual_sum = 0;
		}
	}
	beg = read_beg > 0 ? read_beg : 0;
	for(i = beg; i < read_beg + read_length; ++i) {
		++cinfo[i].depth;
		cinfo[i].mqual_sum += (uint32_t)qual;
	}
	return cinfo;
}

int8_t check_var (faidx_t* fai,
				bamFile fp,
				bam_index_t* idx,
				bam_header_t* header,
				reads* r,
				p_info* cinfo,
				int32_t tid,
				int32_t region_begin,
				int32_t region_end,
				int32_t window_begin,
				int32_t window_end,
				int32_t size,
				int32_t count) {

	int8_t candidate = 0;
	int32_t ref_len;
	char* ref_seq = faidx_fetch_seq(fai, header->target_name[tid], window_begin, window_end, &ref_len);

	if (ref_seq == 0 || ref_len < 1) {
		fprintf(stderr, "Retrieval of reference region \"%s:%d-%d\" failed due to truncated file or corrupt reference index file\n", header->target_name[tid], window_begin, window_end);
		return 0;
	}

	if (pileup_check(fp, idx, ref_seq, tid, window_begin, window_end) > 2) {
		r->count = count;
		call_var (header, r, cinfo, ref_seq, tid, window_begin, window_end, region_begin, region_end, size);
		candidate = 1;
	}

	free(ref_seq);
	return candidate;
}

void slide_window_region (faidx_t* fai, 
						  bamFile fp, 
						  bam1_t* bam, 
						  bam_index_t* idx, 
					      bam_header_t* header,
						  int32_t tid, 
						  int32_t region_begin,	// user required region 
						  int32_t region_end,	// user required region 
						  int32_t size) {

	int32_t n = 128, l = 65536, d = 1024, half_len = 0, count = 0, window_begin = -1, window_end = -1;//, small = 1;
	p_info* cinfo = calloc(d, sizeof(p_info));
	reads* r = calloc(1, sizeof(reads));
	r->pos = malloc(n * sizeof(int32_t));
	r->seq_l = malloc(n * sizeof(int32_t));
	r->qual = malloc(n * sizeof(uint8_t));
	r->seqs = malloc(l * sizeof(uint8_t));	// read sequences stored one after another

	// Buffer the reads.
	bam_iter_t bam_iter = bam_iter_query(idx, tid, region_begin, region_end);	
	while (bam_iter_read (fp, bam_iter, bam) > 0) {
		// Record read information.	
		int32_t read_len = bam->core.l_qseq;
		int32_t char_len = read_len/2;

		if (window_begin == -1) {
			window_begin = bam->core.pos > size ? (bam->core.pos - size) : 0;
			if (window_begin < window_end) window_begin = window_end - WINDOW_EDGE*2;
		}

		if (bam->core.pos - window_begin >= WINDOW_SIZE) {
			if(window_end > window_begin && 2*half_len/(window_end - window_begin) >= 5) {	// average read depth > 5
				cinfo = add_depth(cinfo, &d, bam->core.pos - window_begin, bam->core.l_qseq, bam->core.qual);
				buffer_read1(bam, r, window_begin, window_end, &count, &half_len);	
				check_var (fai, fp, idx, header, r, cinfo, tid, region_begin, region_end, window_begin, window_end, size, count);	
		//		r->count = count;

		//		call_var (header, fai, r, cinfo, tid, window_begin, window_end, region_begin, region_end, size);
			}
			free(r->seqs);
			free(r->qual);
			free(r->seq_l);
			free(r->pos);
			free(r);
			free(cinfo);

			n = 128, l = 65536, d = 1024, half_len = 0, count = 0;
			cinfo = calloc(d, sizeof(p_info));
			r = calloc(1, sizeof(reads));
			r->pos = malloc(n * sizeof(int32_t));
			r->seq_l = malloc(n * sizeof(int32_t));
			r->qual = malloc(n * sizeof(uint8_t));
			r->seqs = malloc(l * sizeof(uint8_t));	// read sequences stored one after another

			window_begin = bam->core.pos > size ? (bam->core.pos - size) : 0;
			if (window_begin < window_end) window_begin = window_end - WINDOW_EDGE*2;
			cinfo = add_depth(cinfo, &d, bam->core.pos - window_begin, bam->core.l_qseq, bam->core.qual);
		} 
		if (bam->core.n_cigar == 0 || bam->core.qual < 10) continue;	// Skip the read that is wrongly mapped or has low mapping quality.
	
		// Adjust memory.
		if (count + 1 >= n) {
			++n;
			kroundup32(n);
			r->pos = realloc(r->pos, n * sizeof(int32_t));	
			r->seq_l = realloc(r->seq_l, n * sizeof(int32_t));	
			r->qual = realloc(r->qual, n * sizeof(uint8_t));
		}
		if (half_len + char_len + 2 >= l) {
			++l;
			kroundup32(l);
			r->seqs = realloc(r->seqs, l * sizeof(uint8_t));
		}
		
		window_end = bam->core.pos + read_len + size;

		cinfo = add_depth(cinfo, &d, bam->core.pos - window_begin, bam->core.l_qseq, bam->core.qual);
		buffer_read1(bam, r, window_begin, window_end, &count, &half_len);
	}

	if(2*half_len/(window_end - window_begin) >= 5) {	// average read depth > 5
				check_var (fai, fp, idx, header, r, cinfo, tid, region_begin, region_end, window_begin, window_end, size, count);	
//		r->count = count;
//		call_var (header, fai, r, cinfo, tid, window_begin, window_end, region_begin, region_end, size);
	}	

	free(r->seqs);
	free(r->qual);
	free(r->seq_l);
	free(r->pos);
	free(r);
	free(cinfo);
}

void slide_window_whole (faidx_t* fai, bamFile fp, bam_header_t* header, bam1_t* bam, bam_index_t* idx, int32_t size) {
	int32_t n = 128, l = 65536, d = 1024, half_len = 0, count = 0, window_begin = -1, window_end = -1, tid = -1;
	p_info* cinfo = calloc(d, sizeof(p_info));
	reads* r = calloc(1, sizeof(reads));
	r->pos = malloc(n * sizeof(int32_t));
	r->seq_l = malloc(n * sizeof(int32_t));
	r->qual = malloc(n * sizeof(uint8_t));
	r->seqs = malloc(l * sizeof(uint8_t));	// read sequences stored one after another

	// Buffer the reads.
	while(bam_read1(fp, bam) > 0){
		// Record read information.	
		int32_t read_len = bam->core.l_qseq;
		int32_t char_len = read_len/2;

		if (window_begin == -1) {
			window_begin = bam->core.pos > size ? (bam->core.pos - size) : 0;
			if (window_begin < window_end) window_begin = window_end - WINDOW_EDGE*2;
			tid = bam->core.tid;
		}

		if ((bam->core.tid != tid) || (bam->core.pos - window_begin >= WINDOW_SIZE)) {
			if(window_end > window_begin && 2*half_len/(window_end - window_begin) >= 5) {	// average read depth > 5
				cinfo = add_depth(cinfo, &d, bam->core.pos - window_begin, bam->core.l_qseq, bam->core.qual);
				buffer_read1(bam, r, window_begin, window_end, &count, &half_len);		
				check_var (fai, fp, idx, header, r, cinfo, tid, -1, 2147483647, window_begin, window_end, size, count);	
			//	r->count = count;
			//	call_var (header, fai, r, cinfo, tid, window_begin, window_end, -1, 2147483647, size);
			}
			free(r->seqs);
			free(r->qual);
			free(r->seq_l);
			free(r->pos);
			free(r);
			free(cinfo);

			n = 128, l = 65536, d = 1024, half_len = 0, count = 0;
			cinfo = calloc(d, sizeof(p_info));
			r = calloc(1, sizeof(reads));
			r->pos = malloc(n * sizeof(int32_t));
			r->seq_l = malloc(n * sizeof(int32_t));
			r->qual = malloc(n * sizeof(uint8_t));
			r->seqs = malloc(l * sizeof(uint8_t));	// read sequences stored one after another

			window_begin = bam->core.pos > size ? (bam->core.pos - size) : 0;
			if ((bam->core.tid == tid) && (window_begin < window_end)) window_begin = window_end - WINDOW_EDGE*2;
			tid = bam->core.tid;
			cinfo = add_depth(cinfo, &d, bam->core.pos - window_begin, bam->core.l_qseq, bam->core.qual);
		} 

		if (bam->core.n_cigar == 0 || bam->core.qual < 10) continue;	// Skip the read that is wrongly mapped or has low mapping quality.
	
		// Adjust memory.
		if (count + 2 >= n) {
			++n;
			kroundup32(n);
			r->pos = realloc(r->pos, n * sizeof(int32_t));	
			r->qual = realloc(r->qual, n * sizeof(uint8_t));
			r->seq_l = realloc(r->seq_l, n * sizeof(int32_t));	
		}
		if (half_len + char_len + 2 >= l) {
			++l;
			kroundup32(l);
			r->seqs = realloc(r->seqs, l * sizeof(uint8_t));
		}
		
		window_end = bam->core.pos + read_len + size;
		
		cinfo = add_depth(cinfo, &d, bam->core.pos - window_begin, bam->core.l_qseq, bam->core.qual);
		buffer_read1(bam, r, window_begin, window_end, &count, &half_len);
	}

	if(2*half_len/(window_end - window_begin) >= 5) {	// average read depth > 5
				check_var (fai, fp, idx, header, r, cinfo, tid, -1, 2147483647, window_begin, window_end, size, count);	
		//r->count = count;
	//	call_var (header, fai, r, cinfo, tid, window_begin, window_end, -1, 2147483647, size);
	}

	free(r->seqs);
	free(r->qual);
	free(r->seq_l);
	free(r->pos);
	free(r);
	free(cinfo);
}

