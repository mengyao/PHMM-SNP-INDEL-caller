/*
 * region.c: Get reference and alignments in a region using samtools-0.1.18
 * Author: Mengyao Zhao
 * Create date: 2011-06-05
 * Last revise data: 2012-12-21
 * Contact: zhangmp@bc.edu 
 */

#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "bam.h"
#include "faidx.h"
#include "hmm.h"
#include "sicall.h"

/*! @function
  @abstract  Round an integer to the next closest power-2 integer.
  @param  x  integer to be rounded (in place)
  @discussion x will be modified.
 */
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

typedef struct {
	double** transition;
	double** emission;
} profile;

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   region [options] <in.fasta> <in.bam> [region1 [...]]\n\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "\t-s N\tN is the largest detectable INDEL length. [default: 100]\n");
	fprintf(stderr, "Notes:\n\
\n\
     A region should be presented in one of the following formats:\n\
     `chr1', `chr2:1,000' and `chr3:1,000-2,000'. When a region is\n\
     specified, the input alignment file must be an indexed BAM file.\n\
\n");
	return 1;
}

// Buffer reads within a 1K window and do Baum-Welch training.
profile* train (int32_t tid,	// reference ID
				char* ref_seq,
				int32_t ref_len,
		   		char* ref_name, 
		   		bamFile* fp,
		   		bam1_t* bam, 
		   		bam_index_t* idx, 
		   		int32_t window_begin, 
		   		int32_t size) {	// maximal detectable INDEL size
	int32_t n = 128, l = 65536, half_len = 0, count = 0; 
	int32_t window_end = window_begin + ref_len + size - 1;
	profile* hmm = (profile*)malloc(sizeof(profile));

	reads* r = calloc(1, sizeof(reads));
	r->pos = malloc(n * sizeof(int32_t));
	r->seq_l = malloc(n * sizeof(int32_t));
	r->seqs = malloc(l * sizeof(uint8_t));	// read sequences stored one after another

	// Retrieve the alignments that are overlapped with the specified region.
	bam_iter_t bam_iter = bam_iter_query(idx, tid, window_begin, window_end);	
	while (bam_iter_read (*fp, bam_iter, bam) >= 0) {

		if (bam->core.n_cigar == 0) continue;	// Skip the read that is wrongly mapped.
		// Record read information.
		uint32_t* cigar = bam1_cigar(bam);
		uint8_t* read_seq = bam1_seq(bam);
		int32_t read_len = bam->core.l_qseq, j;
		int32_t char_len = read_len/2;

		// Adjust memory.
		if (count + 1 >= n) {
			++n;
			kroundup32(n);
			r->pos = realloc(r->pos, n * sizeof(int32_t));	
			r->seq_l = realloc(r->seq_l, n * sizeof(int32_t));	
		}
		if (half_len + char_len + 1 >= l) {
			++l;
			kroundup32(l);
			r->seqs = realloc(r->seqs, l * sizeof(uint8_t));
		}
	
		// Buffer reads in this region.
		if (bam->core.pos < window_begin) {	
		// read head is aligned out of the window: truncate the head
			int32_t pos = bam->core.pos;
			int32_t clip_len = 0;
			while (pos < window_begin) {
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
			}
			if (clip_len%2) {	// Remove one more read sequence residual.
				read_len -= (clip_len + 1);
				if (read_len == 0) continue;
				read_seq += ((clip_len + 1)/2);
				r->pos[count] = window_begin + 1;
			} else {
				read_len -= clip_len;
				read_seq += (clip_len/2);
				r->pos[count] = window_begin;
			}
		} else if (bam->core.pos + read_len > window_end) {	
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
		}	
		if (bam->core.pos >= window_begin) r->pos[count] = bam->core.pos;
		r->seq_l[count] = read_len;
		char_len = read_len/2;
		for (j = half_len; j < half_len + char_len; j ++) r->seqs[j] = read_seq[j - half_len];
		if (read_len%2) r->seqs[j] = read_seq[j - half_len];
		half_len += char_len;
		if (read_len%2) half_len ++;
		count ++;
	}

	if (count == 0) hmm = NULL;
	else {
		hmm->transition = transition_init (0.002, 0.98, 0.00067, 0.02, 0.998, ref_len + size);
		hmm->emission = emission_init(ref_seq, size);
		r->count = count;
		baum_welch (hmm->transition, hmm->emission, ref_seq, window_begin, ref_len + size, size, r, 0.01); /* 0-based coordinate */ 
	}

	bam_iter_destroy(bam_iter);
	free(r->seqs);
	free(r->seq_l);
	free(r->pos);
	free(r);
	
	return hmm;
}

// Deal with one large region with sliding window.
int32_t slide_window (faidx_t* fai, 
					  bam_header_t* header, 
					  bamFile fp, 
					  bam1_t* bam, 
					  bam_index_t* idx, 
					  int32_t tid, 
					  int32_t beg,
					  int32_t end,	// When end = 0, it means the slide_window will deal with the whole chromosome. 
					  int32_t size) {

	int32_t ref_len, window_end; 
	char* ref_seq;
	profile* hmm;
	if (end == 0 || end - beg > 1000) {	// large region or whole genome (genome length >= 1000)
		window_end = end == 0 ? beg + 999 : (beg + 999 > end ? end : beg + 999);
		ref_seq = faidx_fetch_seq(fai, header->target_name[tid], beg, window_end, &ref_len);
		if (ref_seq == 0 || ref_len < 1) {
			fprintf(stderr, "Retrieval of reference region \"%s:%d-%d\" failed due to truncated file or corrupt reference index file\n", header->target_name[tid], beg, beg + 999);
			return 0;
		}
		while (ref_seq != 0 && ref_len >= 100) {
			hmm = train(tid, ref_seq, ref_len, header->target_name[tid], &fp, bam, idx, beg, size);
			if (hmm) {
				int32_t region_beg = beg + 50;
				int32_t region_end = beg + ref_len - 50;
				likelihood (hmm->transition, hmm->emission, ref_seq, header->target_name[tid], beg, region_beg, region_end, 0);	
				transition_destroy(hmm->transition, ref_len + size);
				emission_destroy(hmm->emission, ref_len + size);
				free(hmm);
			}
			free(ref_seq);
			beg += 900;	
			window_end = end == 0 ? beg + 999 : (beg + 999 > end ? end : beg + 999);
			if (window_end > beg) ref_seq = faidx_fetch_seq(fai, header->target_name[tid], beg, window_end, &ref_len);
			else break;
		}
	} else {	// small region
		ref_seq = faidx_fetch_seq(fai, header->target_name[tid], beg, end, &ref_len);
		if (ref_seq == 0 || ref_len < 1) {
			fprintf(stderr, "Retrieval of reference region \"%s:%d-%d\" failed due to truncated file or corrupt reference index file.\n", header->target_name[tid], beg, end);
			return 0;
		}
		hmm = train(tid, ref_seq, ref_len, header->target_name[tid], &fp, bam, idx, beg, size);
		if (hmm) {
			int32_t region_beg = beg + ref_len / 10;
			int32_t region_end = beg + beg + 9 * ref_len / 10 - 1;
			likelihood (hmm->transition, hmm->emission, ref_seq, header->target_name[tid], beg, region_beg, region_end, 0);	
			transition_destroy(hmm->transition, ref_len + size);
			emission_destroy(hmm->emission, ref_len + size);
			free(hmm);
		}
		else fprintf(stderr, "There's no read falling into the given region \"%s:%d-%d\".\n", header->target_name[tid], beg, end);	
		free(ref_seq);
	}	
	return 1;
}

int main (int argc, char * const argv[]) {
	int32_t ret = 0;	// return value
	float cpu_time;
	clock_t start, end;
	start = clock();

	int32_t l,i;
	int32_t size = 101;	// default largest detectable INDEL length
	bamFile fp;
	bam_header_t* header;
	bam1_t* bam = bam_init1();
	bam_index_t* idx = 0;
	faidx_t* fai;

	fprintf (stdout, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

	// Parse command line.
	while ((l = getopt(argc, argv, "s:")) >= 0) {
		switch (l) {
			case 's': size = atoi(optarg); break;
		}
	}

	if (optind + 2 > argc) return usage();

	if ((fai = fai_load(argv[optind])) == 0) {
		fprintf(stderr, "Random alignment retrieval requires the reference index file.\n");
		ret = 1;
		goto end_fai;
	}
	if ((fp = bam_open(argv[optind + 1], "r")) == 0) {
		fprintf(stderr, "Fail to open \"%s\" for reading.\n", argv[2]);
		ret = 1;
		goto end_fp;
	}
	if ((header = bam_header_read(fp)) == 0) {
		fprintf(stderr, "Fail to read the header from \"%s\".\n", argv[2]);
		ret = 1;
		goto end_header;
	}
	if ((idx = bam_index_load(argv[optind + 1])) == 0) { // index is unavailable
		fprintf(stderr, "Random alignment retrieval only works for indexed BAM files.\n");
		ret = 1;
		goto end_idx;
	}

	if (argc == (optind + 2)) {	// No region is given by the command line.
		int32_t ref_count = faidx_fetch_nseq(fai), tid;
		for (tid = 0; tid < ref_count; ++tid) if (! slide_window (fai, header, fp, bam , idx, tid, 0, 0, size)) continue; 
	} else {	// Regions are given by the command line.
		int32_t tid, beg, end;
		i = optind + 2;
		while(i < argc) {
			bam_parse_region(header, argv[i], &tid, &beg, &end); // parse a region in the format like `chr2:100-200'
			if (tid < 0) { // reference name is not found
				fprintf(stderr, "region \"%s\" specifies an unknown reference name.\n", argv[i]);
				ret = 1;
				goto end;
			}
			if (! slide_window (fai, header, fp, bam , idx, tid, beg, end, size)) continue; 
			++i;
		}
	}
end:
	bam_index_destroy(idx);
end_idx:
	bam_header_destroy(header);
end_header:
	bam_close(fp);
end_fp:
	fai_destroy(fai);
end_fai:
	bam_destroy1(bam);

	end = clock();
	cpu_time = ((float) (end - start)) / CLOCKS_PER_SEC;
	fprintf(stdout, "\n\nCPU time: %f seconds\n", cpu_time);

	return ret;
}

