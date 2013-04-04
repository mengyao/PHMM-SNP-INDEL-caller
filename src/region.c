/*
 * region.c: Get reference and alignments in a region using samtools-0.1.18
 * Author: Mengyao Zhao
 * Create date: 2011-06-05
 * Last revise date: 2013-04-04
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
#define WINDOW_EDGE 50

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
		if (read_len < 13) return 0;
	}	
	if (bam->core.pos >= window_begin) r->pos[*count] = bam->core.pos;
	r->seq_l[*count] = read_len;
	char_len = read_len/2;
	for (j = *half_len; j < *half_len + char_len; j ++) r->seqs[j] = read_seq[j - *half_len];
	if (read_len%2) r->seqs[j] = read_seq[j - *half_len];
	(*half_len) += char_len;
	if (read_len%2) (*half_len) ++;
	(*count) ++;
	
	return 1;
}

// Buffer reads within a 1K window and do Baum-Welch training.
void train (int32_t tid,	// reference ID
				char* ref_seq,
				int32_t ref_len,
		   		char* ref_name,
				bam_header_t* header, 
		   		bamFile fp,
		   		bam1_t* bam, 
		   		bam_index_t* idx, 
		   		int32_t window_begin, 
		   		int32_t size) {	// maximal detectable INDEL size

	int32_t n = 128, l = 65536, d = 1024, half_len = 0, count = 0, i, beg; 
	int32_t window_end = window_begin + ref_len + size - 1;
	profile* hmm = (profile*)malloc(sizeof(profile));
	int32_t* depth = calloc(d, sizeof(int32_t));

	reads* r = calloc(1, sizeof(reads));
	r->pos = malloc(n * sizeof(int32_t));
	r->seq_l = malloc(n * sizeof(int32_t));
	r->seqs = malloc(l * sizeof(uint8_t));	// read sequences stored one after another

	// Retrieve the alignments that are overlapped with the specified region.
	bam_iter_t bam_iter = bam_iter_query(idx, tid, window_begin, window_end);	
	while (bam_iter_read (fp, bam_iter, bam) >= 0) {
		if (bam->core.n_cigar == 0) continue;	// Skip the read that is wrongly mapped.
		// Record read information.
		int32_t read_len = bam->core.l_qseq;
		int32_t char_len = read_len/2;

		// Adjust memory.
		if(bam->core.pos + read_len - window_begin > d) {
			++d;
			kroundup32(d);
			depth = realloc(depth, d*sizeof(int32_t));
		}
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

		// count and half_len are updated in this function.
		if(buffer_read1(bam, r, window_begin, window_end, &count, &half_len)) {
			beg = bam->core.pos - window_begin > 0 ? bam->core.pos - window_begin : 0;
			for(i = beg; i < bam->core.pos + bam->core.l_qseq - window_begin; ++i) {
			//	fprintf(stderr, "i: %d\n", i);
				++depth[i];
			}
		}else continue;		
	}

	if (count == 0) fprintf(stderr, "There's no read falling into the given region \"%s:%d-%d\".\n", header->target_name[tid], window_begin, window_end);	
	else {
		int32_t frame_begin = window_begin + ref_len / 10;
		int32_t frame_end = window_begin + 9 * ref_len / 10 - 1;
		hmm->transition = transition_init (0.002, 0.98, 0.00067, 0.02, 0.998, ref_len + size);
		hmm->emission = emission_init(ref_seq, size);
		r->count = count;
		baum_welch (hmm->transition, hmm->emission, ref_seq, window_begin, ref_len + size, size, r, 0.01); /* 0-based coordinate */
		likelihood (header, hmm->transition, hmm->emission, ref_seq, depth, tid, window_begin, frame_begin, frame_end, size, 0);
		transition_destroy(hmm->transition, ref_len + size);
		emission_destroy(hmm->emission, ref_len + size);
		free(hmm);
	}

	bam_iter_destroy(bam_iter);
	free(r->seqs);
	free(r->seq_l);
	free(r->pos);
	free(r);
	free(depth);
}

void call_var (//bamFile fp,
			   bam_header_t* header,
			 //  bam_index_t* idx,
				faidx_t* fai,
				  reads* r, 
					int32_t* depth,
			   	  int32_t tid, 
			   	  int32_t window_begin, 
			   	  int32_t window_end,
			   	  int32_t region_begin,	// only used in slide_window_region
			   	  int32_t region_end,	// only used in slide_window_region 
			   	  int32_t size) {

	int32_t ref_len, frame_begin, frame_end, temp;
	char* ref_seq = faidx_fetch_seq(fai, header->target_name[tid], window_begin, window_end, &ref_len);
	profile* hmm = (profile*)malloc(sizeof(profile));

	if (ref_seq == 0 || ref_len < 1) {
		fprintf(stderr, "Retrieval of reference region \"%s:%d-%d\" failed due to truncated file or corrupt reference index file\n", header->target_name[tid], window_begin, window_end);
		return;
	}

	hmm->transition = transition_init (0.002, 0.98, 0.00067, 0.02, 0.998, ref_len + size);
	hmm->emission = emission_init(ref_seq, size);
	baum_welch (hmm->transition, hmm->emission, ref_seq, window_begin, ref_len + size, size, r, 0.01); 
 
	temp = window_begin + WINDOW_EDGE;
	frame_begin = temp > region_begin ? temp : region_begin;
	temp = window_begin + ref_len - WINDOW_EDGE;
	frame_end = temp < region_end ? temp : region_end;
	if(frame_end > frame_begin) 
		likelihood (header, hmm->transition, hmm->emission, ref_seq, depth, tid, window_begin, frame_begin, frame_end, size, 0);	
	transition_destroy(hmm->transition, ref_len + size);
	emission_destroy(hmm->emission, ref_len + size);
	free(hmm);
	free(ref_seq);

	return;
}

void slide_window_region (faidx_t* fai, 
						  bamFile fp, 
						  bam1_t* bam, 
						  bam_index_t* idx, 
					      bam_header_t* header,
						  int32_t tid, 
						  int32_t region_begin, 
						  int32_t region_end, 
						  int32_t size) {

	int32_t n = 128, l = 65536, d = 1024, half_len = 0, count = 0, window_begin = -1, window_end = -1, i, beg;
	int32_t* depth = calloc(d, sizeof(int32_t));
	reads* r = calloc(1, sizeof(reads));
	r->pos = malloc(n * sizeof(int32_t));
	r->seq_l = malloc(n * sizeof(int32_t));
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

		if ((bam->core.pos - window_begin >= 1000) && (count >= 100)) {
			if(2*half_len/(window_end - window_begin - 2*size) > 5) {	// average read depth > 5
				beg = bam->core.pos - window_begin > 0 ? bam->core.pos - window_begin : 0;
				for(i = beg; i < bam->core.pos + bam->core.l_qseq - window_begin; ++i) ++depth[i];
				buffer_read1(bam, r, window_begin, window_end, &count, &half_len);		
				r->count = count;
				call_var (header, fai, r, depth, tid, window_begin, window_end, region_begin, region_end, size);
			}
			free(r->seqs);
			free(r->seq_l);
			free(r->pos);
			free(r);
			free(depth);

			n = 128, l = 65536, d = 1024, half_len = 0, count = 0;
			depth = calloc(d, sizeof(int32_t));
			r = calloc(1, sizeof(reads));
			r->pos = malloc(n * sizeof(int32_t));
			r->seq_l = malloc(n * sizeof(int32_t));
			r->seqs = malloc(l * sizeof(uint8_t));	// read sequences stored one after another

			window_begin = bam->core.pos > size ? (bam->core.pos - size) : 0;
			if (window_begin < window_end) window_begin = window_end - WINDOW_EDGE*2;

			beg = bam->core.pos - window_begin > 0 ? bam->core.pos - window_begin : 0;
			for(i = beg; i < bam->core.pos + bam->core.l_qseq - window_begin; ++i) ++depth[i];
			
			// Buffer the information of one read.
			buffer_read1(bam, r, window_begin, window_end, &count, &half_len);		
		} 

		if (bam->core.n_cigar == 0) continue;	// Skip the read that is wrongly mapped.
	
		// Adjust memory.
		if(bam->core.pos + read_len - window_begin + 1 > d) {
			++d;
			kroundup32(d);
			depth = realloc(depth, d*sizeof(int32_t));
		}
		if (count + 1 >= n) {
			++n;
			kroundup32(n);
			r->pos = realloc(r->pos, n * sizeof(int32_t));	
			r->seq_l = realloc(r->seq_l, n * sizeof(int32_t));	
		}
		if (half_len + char_len + 2 >= l) {
			++l;
			kroundup32(l);
			r->seqs = realloc(r->seqs, l * sizeof(uint8_t));
		}
		
		window_end = bam->core.pos + read_len + size;

		beg = bam->core.pos - window_begin > 0 ? bam->core.pos - window_begin : 0;
		for(i = beg; i < bam->core.pos + bam->core.l_qseq - window_begin; ++i) ++depth[i];
		// Buffer the information of one read. Skip, if the read length turns to 0 after truncation.
		buffer_read1(bam, r, window_begin, window_end, &count, &half_len);
	}

	if(2*half_len/(window_end - window_begin - 2*size) > 5) {	// average read depth > 5
		r->count = count;
		call_var (header, fai, r, depth, tid, window_begin, window_end, region_begin, region_end, size);
	}

	free(r->seqs);
	free(r->seq_l);
	free(r->pos);
	free(r);
	free(depth);
}

// Deal with one user request region.
int32_t region(faidx_t* fai, 
					  bam_header_t* header, 
					  bamFile fp, 
					  bam1_t* bam, 
					  bam_index_t* idx, 
					  char* region_str,
					  int32_t size) {

	int32_t tid, ref_len, region_begin, region_end;
	char* ref_seq;

	bam_parse_region(header, region_str, &tid, &region_begin, &region_end); // parse a region in the format like `chr2:100-200'
	if (tid < 0) { // reference name is not found
		fprintf(stderr, "region \"%s\" specifies an unknown reference name.\n", region_str);
		return 0;
	}

	if (region_end - region_begin > 1000) 	// large region (genome length >= 1000)
		slide_window_region(fai, fp, bam, idx, header, tid, region_begin, region_end, size);
	else {	// small region
		ref_seq = faidx_fetch_seq(fai, header->target_name[tid], region_begin, region_end, &ref_len);
		if (ref_seq == 0 || ref_len < 1) {
			fprintf(stderr, "Retrieval of reference region \"%s:%d-%d\" failed due to truncated file or corrupt reference index file.\n", header->target_name[tid], region_begin, region_end);
			return 0;
		}
		train(tid, ref_seq, ref_len, header->target_name[tid], header, fp, bam, idx, region_begin, size);
		free(ref_seq);
	}	
	return 1;
}

void slide_window_whole (faidx_t* fai, bamFile fp, bam_header_t* header, bam1_t* bam, bam_index_t* idx, int32_t size) {
//	int32_t window_end = -1, tid = -1, one_read = 0, bam_end = 1;
 
//	while (bam_end > 0) {
	int32_t n = 128, l = 65536, d = 1024, half_len = 0, count = 0, window_begin = -1, window_end = -1, tid = -1, i, beg;
	int32_t* depth = calloc(d, sizeof(int32_t));
	reads* r = calloc(1, sizeof(reads));
	r->pos = malloc(n * sizeof(int32_t));
	r->seq_l = malloc(n * sizeof(int32_t));
	r->seqs = malloc(l * sizeof(uint8_t));	// read sequences stored one after another

//	if (one_read == 1) {	// the 1st read in the new window		
//	}

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

		if ((bam->core.tid != tid) || ((bam->core.pos - window_begin >= 1000) && (count >= 100))) {
			if(2*half_len/(window_end - window_begin - 2*size) > 5) {	// average read depth > 5
				beg = bam->core.pos - window_begin > 0 ? bam->core.pos - window_begin : 0;
				for(i = beg; i < bam->core.pos + bam->core.l_qseq - window_begin; ++i) ++depth[i];
				buffer_read1(bam, r, window_begin, window_end, &count, &half_len);		
				r->count = count;
				call_var (header, fai, r, depth, tid, window_begin, window_end, -1, 2147483647, size);
			}
			free(r->seqs);
			free(r->seq_l);
			free(r->pos);
			free(r);
			free(depth);
//			one_read = 1;
//			break;	// Close a window.
			n = 128, l = 65536, d = 1024, half_len = 0, count = 0;
			depth = calloc(d, sizeof(int32_t));
			r = calloc(1, sizeof(reads));
			r->pos = malloc(n * sizeof(int32_t));
			r->seq_l = malloc(n * sizeof(int32_t));
			r->seqs = malloc(l * sizeof(uint8_t));	// read sequences stored one after another

			window_begin = bam->core.pos > size ? (bam->core.pos - size) : 0;
			if ((bam->core.tid == tid) && (window_begin < window_end)) window_begin = window_end - WINDOW_EDGE*2;
			tid = bam->core.tid;

			beg = bam->core.pos - window_begin > 0 ? bam->core.pos - window_begin : 0;
			for(i = beg; i < bam->core.pos + bam->core.l_qseq - window_begin; ++i) ++depth[i];
			buffer_read1(bam, r, window_begin, window_end, &count, &half_len);		
		} 

		if (bam->core.n_cigar == 0) continue;	// Skip the read that is wrongly mapped.
	
		// Adjust memory.
		if(bam->core.pos + read_len - window_begin > d) {
			++d;
			kroundup32(d);
			depth = realloc(depth, d*sizeof(int32_t));
		}
		if (count + 2 >= n) {
			++n;
			kroundup32(n);
			r->pos = realloc(r->pos, n * sizeof(int32_t));	
			r->seq_l = realloc(r->seq_l, n * sizeof(int32_t));	
		}
		if (half_len + char_len + 2 >= l) {
			++l;
			kroundup32(l);
			r->seqs = realloc(r->seqs, l * sizeof(uint8_t));
		}
		
		window_end = bam->core.pos + read_len + size;
		
		// Buffer the information of one read. Skip, if the read length turns to 0 after truncation.
		beg = bam->core.pos - window_begin > 0 ? bam->core.pos - window_begin : 0;
		for(i = beg; i < bam->core.pos + bam->core.l_qseq - window_begin; ++i) ++depth[i];
		buffer_read1(bam, r, window_begin, window_end, &count, &half_len);
	}

	if(2*half_len/(window_end - window_begin - 2*size) > 5) {	// average read depth > 5
//		beg = bam->core.pos - window_begin > 0 ? bam->core.pos - window_begin : 0;
//		for(i = beg; i < bam->core.pos + bam->core.l_qseq - window_begin; ++i) ++depth[i];
//		buffer_read1(bam, r, window_begin, window_end, &count, &half_len);		
		r->count = count;
		call_var (header, fai, r, depth, tid, window_begin, window_end, -1, 2147483647, size);
	}

	free(r->seqs);
	free(r->seq_l);
	free(r->pos);
	free(r);
	free(depth);
//	}
}

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

	if (argc == (optind + 2)) slide_window_whole(fai, fp, header, bam, idx, size); 	// No region is given by the command line.
	else {	// Regions are given by the command line.
		i = optind + 2;
		while(i < argc) {
			region(fai, header, fp, bam, idx, argv[i], size); 
			++i;
		}
	}

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

