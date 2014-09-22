/*
 * pileup.c: Evaluate the regions with candidate variations from the pileup data using samtools-0.1.18
 * Author: Mengyao Zhao
 * Create date: 2014-07-23
 * Last revise date: 2014-08-11
 * Contact: zhangmp@bc.edu 
 */

#include <stdio.h>
#include "pileup.h"
#include "viterbi.h"

typedef struct {     // auxiliary data structure
	bamFile fp;      // the file handler
	bam_iter_t iter; // NULL if a region not specified
	int min_mapQ;    // mapQ filter
} aux_t;

// This function reads a BAM alignment from one BAM file.
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
	aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
	int ret = aux->iter? bam_iter_read(aux->fp, aux->iter, b) : bam_read1(aux->fp, b);
	if ((int)b->core.qual < aux->min_mapQ) b->core.flag |= BAM_FUNMAP;
	return ret;
}

int32_t pileup_check (bamFile fp, 
			 		   bam_index_t* idx,
						char* ref_seq, 
			           int32_t tid, 
			  	       int32_t beg,
				   	int32_t end) {

//	fprintf(stderr, "beg: %d\tend: %d\n", beg, end);
	// Get the average base pileup around the candidate variation location.
	int* n_plp = calloc(1, sizeof(int)); 	// the number of covering reads from the BAM file
	int32_t  mapQ = 20, baseQ = 9, pos = beg, max = 0;
	const bam_pileup1_t **plp = calloc(1, sizeof(void*)); 	// points to the array of covering reads (internal in mplp)
	aux_t *data;
	bam_mplp_t mplp;
	
	// initialize the auxiliary data structures
	data = calloc(1, sizeof(aux_t));
	data->fp = fp;
	data->min_mapQ = mapQ;                    // set the mapQ filter
	data->iter = bam_iter_query(idx, tid, beg, end); // set the iterator


	//fprintf(stderr, "%p,%p,%p\n", data, &data, &mplp); data[0]
	// the core multi-pileup loop
	mplp = bam_mplp_init(1, read_bam, (void**)&data); // initialization
	while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) { // come to the next covered position
		if (pos < beg || pos >= end) continue; // out of range; skip
//	fprintf(stderr, "tid: %d\tpos: %d\tbeg: %d\tend: %d\n", tid, pos, beg, end);
		int j, m = 0;
		for (j = 0; j < n_plp[0]; ++j) {
			const bam_pileup1_t *p = plp[0] + j;
			const uint8_t *seq = bam1_seq(p->b); 
			if (bam1_qual(p->b)[p->qpos] > baseQ && p->indel) { 
				++m; // having indels
			//	fprintf(stderr, "m: %d\n", m);
			}
		//	printf("%d\n", bam1_seqi(seq, p->qpos));
			if (bam1_qual(p->b)[p->qpos] > baseQ && (! p->is_del) && num2base[bam1_seqi(seq, p->qpos)] != ref_seq[pos - beg]) {
				++m;// printf("%d\n", bam1_seqi(seq, p->qpos)); // print each base
//				fprintf(stderr, "read: %c\tref: %c\n", num2base[bam1_seqi(seq, p->qpos)], ref_seq[pos - beg]);
			}
		}
//		if (m >= 2) fprintf(stderr, "*pos: %d*\tm: %d\n", pos, m);
		if (m > max) max = m;
	}
	free(n_plp); free(plp);
	bam_mplp_destroy(mplp);
	return max;
}

