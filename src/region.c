/*
 * region.c: Get reference and alignments in a region using samtools-0.1.18
 * Author: Mengyao Zhao
 * Create date: 2011-06-05
 * Last revise data: 2012-10-24
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

		int32_t n = 128, l = 65536, c = 4096, half_len = 0, count = 0, cigar_len = 0, window_end = window_begin + 999;
		profile* hmm = (profile*)malloc(sizeof(profile));;

		reads* r = calloc(1, sizeof(reads));
		r->pos = malloc(n * sizeof(int32_t));
		r->seq_l = malloc(n * sizeof(int32_t));
		r->cigar = malloc(c * sizeof(uint32_t));	// cigar strings of reads stored one after another
      //  memset(r->cigar, 0, c * sizeof (uint32_t));
		r->n_cigar = malloc(n * sizeof(uint16_t));	// length of cigar string
		r->seqs = malloc(l * sizeof(uint8_t));	// read sequences stored one after another
	
		// Retrieve the alignments that are overlapped with the specified region.
//		fprintf(stderr, "beg: %d\nend: %d\ntid: %d\n", beg, end, tid);	
		bam_iter_t bam_iter = bam_iter_query(idx, tid, window_begin, window_end);	
		
//		fprintf(stderr, "bam_iter: %d\nfp: %d\nbam: %d\n", bam_iter, fp, bam);

//FIXME: move read truncation into this loop
		while (bam_iter_read (*fp, bam_iter, bam) >= 0) {
//				fprintf(stderr, "here!!!\n");
			uint32_t* cigar = bam1_cigar(bam);
			uint8_t* read_seq = bam1_seq(bam);
			int32_t read_len = bam->core.l_qseq;
			int32_t char_len = read_len/2, j;
		//	have_read = 1;
			if (count + 1 >= n) {
				++n;
				kroundup32(n);
				r->pos = realloc(r->pos, n * sizeof(int32_t));	
				r->seq_l = realloc(r->seq_l, n * sizeof(int32_t));	
				r->n_cigar = realloc(r->n_cigar, n * sizeof(int32_t));
			}
			r->pos[count] = bam->core.pos;
			r->seq_l[count] = read_len;
			r->n_cigar[count] = bam->core.n_cigar;
			count ++;

			if (half_len * 2 + char_len + 1 >= l) {
				++l;
				kroundup32(l);
				r->seqs = realloc(r->seqs, l * sizeof(uint8_t));
			}		
			for (j = half_len; j < half_len + char_len; j ++) {
			//	fprintf(stderr, "read_seq[%d]: %d\t", j - half_len, read_seq[j - half_len] );
				r->seqs[j] = read_seq[j - half_len];
			}
			if (read_len%2) r->seqs[j] = read_seq[j - half_len];
			half_len += char_len + read_len%2;
		//	fprintf(stderr, "\n");


			if (cigar_len + bam->core.n_cigar >= c) {
				++c;
				kroundup32(c);
				r->cigar = realloc(r->cigar, c * sizeof(uint32_t));
			}
			for (j = cigar_len; j < cigar_len + bam->core.n_cigar; ++j) r->cigar[j] = cigar[j - cigar_len];
		//	fprintf(stderr, "cigar_len: %d\n", cigar_len);
			cigar_len += bam->core.n_cigar;
		}

		if (count == 0) {
			//fprintf (stderr, "There is no read in the given region \"%s:%d-%d\".\n", ref_name, beg, end);
			//free(hmm);
			hmm = NULL;
		} else {
			hmm->transition = transition_init (0.002, 0.98, 0.00067, 0.02, 0.998, ref_len + size);
			hmm->emission = emission_init(ref_seq, size);
			r->count = count;
			baum_welch (hmm->transition, hmm->emission, ref_seq, window_begin, ref_len + size, size, r, 0.01); /* 0-based coordinate */ 
		}

		bam_iter_destroy(bam_iter);
		free(r->seqs);
		free(r->n_cigar);
		free(r->cigar);
		free(r->seq_l);
		free(r->pos);
		free(r);
		
		return hmm;
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
	bam_index_t *idx = 0;
	faidx_t* fai;
	profile* hmm;

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

	if (argc == (optind + 2)) {
		int ref_count = faidx_fetch_nseq(fai);
		int32_t tid, ref_len;
		char* ref_seq;
		for (tid = 0; tid < ref_count; ++tid) {
			int32_t beg = 0; 
			ref_seq = faidx_fetch_seq(fai, header->target_name[tid], beg, beg + 999, &ref_len);
			if (ref_seq == 0 || ref_len < 1) {
				fprintf(stderr, "Retrieval of reference region \"%s:%d-%d\" failed due to truncated file or corrupt reference index file\n", header->target_name[tid], beg, beg + 999);
			//	ret = 1;
			//	goto end;
				continue;
			}
			while (ref_seq != 0 && ref_len > 0) {
				hmm = train(tid, ref_seq, ref_len, header->target_name[tid], &fp, bam, idx, beg, size);
				if (hmm) {
				//	ret = 1;
				//	goto end;
			//	} else {
					int32_t region_beg = ref_len > 1000 ? beg + 100 : beg + ref_len / 10;
					int32_t region_end = ref_len > 1000 ? beg + 899 : beg + 9 * ref_len / 10 - 1;	
					likelihood (hmm->transition, hmm->emission, ref_seq, header->target_name[tid], beg, region_beg, region_end, 0);	
			//		fprintf(stderr, "ref_name: %s\nbeg: %d\nref_len: %d\nref_seq: %s\n", header->target_name[tid], beg, ref_len, ref_seq);

					transition_destroy(hmm->transition, ref_len + size);
					emission_destroy(hmm->emission, ref_len + size);
				}
				free(hmm);
				free(ref_seq);
				beg += 800;
				ref_seq = faidx_fetch_seq(fai, header->target_name[tid], beg, beg + 999, &ref_len);
			}
			free(ref_seq);
		}
	} else {
	//	fprintf(stderr, "here\n");
	//	fprintf(stderr, "argc: %d\toptind: %d\n", argc, optind);

		int32_t tid, ref_len, beg, end;
		i = optind + 2;
	//	fprintf(stderr, "optind: %d\nargv[i]: %s\n", optind, argv[i]);
		bam_parse_region(header, argv[i], &tid, &beg, &end); // parse a region in the format like `chr2:100-200'
		if (tid < 0) { // reference name is not found
			fprintf(stderr, "region \"%s\" specifies an unknown reference name.\n", argv[i]);
			ret = 1;
			goto end;
		}
		while(i < argc) {
			char* ref_seq;
			beg = beg - 100 > 0 ? beg - 100 : 0;
			ref_seq = faidx_fetch_seq(fai, header->target_name[tid], beg, beg + 999, &ref_len);
			if (ref_seq == 0 || ref_len < 1) {
				fprintf(stderr, "Retrieval of reference region \"%s:%d-%d\" failed due to truncated file or corrupt reference index file\n", header->target_name[tid], beg, beg + 999);
	//			ret = 1;
	//			goto end;
				continue;
			}
			if (ref_seq != 0) {
			hmm = train(tid, ref_seq, ref_len, header->target_name[tid], &fp, bam, idx, beg, size);
				if (! hmm)// {} 
				//	ret = 1;
				//	goto end;
					fprintf (stderr, "There is no read in the given region \"%s:%d-%d\".\n", header->target_name[tid], beg, end);
				else {
					while (end < beg + 1000) {
						likelihood (hmm->transition, hmm->emission, ref_seq, header->target_name[tid], beg, beg + 100, end, 0);	
						if (++i >= argc) break;
						bam_parse_region(header, argv[i], &tid, &beg, &end); // parse a region in the format like `chr2:100-200'
					}
					transition_destroy(hmm->transition, ref_len + size);
					emission_destroy(hmm->emission, ref_len + size);
				}
				free(hmm);
			}
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

