/*
 * sicall.c: Calculate the probabilities of SNPs and INDELs. Make the calling.
 * Author: Mengyao Zhao
 * Create date: 2011-08-09
 * Contact: zhangmp@bc.edu
 * Last revise: 2013-03-13 
 */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "sicall.h"
#include "bam.h"

typedef struct {
	int8_t num;
	double prob;
} p_max;

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

p_max* refp (double** emission, char* ref, int32_t k) {
	p_max* ref_allele = (p_max*)malloc(sizeof(p_max));
	switch (ref[k]) {
		case 'A':
		case 'a':
			ref_allele->prob = emission[k + 1][1];
			ref_allele->num = 1;
			break;
		case 'C':
		case 'c':
			ref_allele->prob = emission[k + 1][2];
			ref_allele->num = 2;
			break;
		case 'G':
		case 'g':
			ref_allele->prob = emission[k + 1][4];
			ref_allele->num = 4;
			break;
		case 'T':
		case 't':
			ref_allele->prob = emission[k + 1][8];
			ref_allele->num = 8;
			break;
		default:
			fprintf(stderr, "Wrong reference sequence. \n");
			exit (1);
			break;
	}
	return ref_allele;
}

p_max* bubble3 (int8_t n1, double p1, int8_t n2, double p2, int8_t n3, double p3) {
	p_max *m = (p_max*)malloc(sizeof(p_max));
	m->num = n1;
	m->prob = p1;
	if (p2 > m->prob) {
		m->num = n2;
		m->prob = p2;
	}
	if (p3 > m->prob) {
		m->num = n3;
		m->prob = p3;
	}
	return m;
}

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
		default:
			fprintf(stderr, "The base number is assigned wrongly.\n");
			exit(1);
	} 
	return base;
}

float base_read_depth (bamFile fp, 
			 		   bam_index_t* idx, 
			           int32_t tid, 
			  	  	   int32_t pos,
			  	       int32_t beg,
				   	int32_t end) {

	// Get the average base pileup around the candidate variation location.
	int8_t i, n = 1;	// There's only one BAM file as input.
	int* n_plp;
	int32_t  mapQ = 0, baseQ = 0, total_depth = 0, ave_depth;
	const bam_pileup1_t **plp;
	aux_t **data;
	bam_mplp_t mplp;
	
	// initialize the auxiliary data structures
	data = calloc(n, sizeof(void*)); // data[i] for the i-th input
	for (i = 0; i < n; ++i) {
		data[i] = calloc(1, sizeof(aux_t));
		data[i]->fp = fp;
		data[i]->min_mapQ = mapQ;                    // set the mapQ filter
		data[i]->iter = bam_iter_query(idx, tid, beg, end); // set the iterator
	}

	// the core multi-pileup loop
	mplp = bam_mplp_init(n, read_bam, (void**)data); // initialization
	n_plp = calloc(n, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
	plp = calloc(n, sizeof(void*)); // plp[i] points to the array of covering reads (internal in mplp)
	while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) { // come to the next covered position
		if (pos < beg || pos >= end) continue; // out of range; skip
		for (i = 0; i < n; ++i) { // base level filters have to go here
			int j, m = 0;
			for (j = 0; j < n_plp[i]; ++j) {
				const bam_pileup1_t *p = plp[i] + j; // DON'T modfity plp[][] unless you really know
				if (p->is_del || p->is_refskip) ++m; // having dels or refskips at tid:pos
				else if (bam1_qual(p->b)[p->qpos] < baseQ) ++m; // low base quality
			}
			total_depth += n_plp[i] - m;	// add up base read depth
		}
	}
	free(n_plp); free(plp);
	bam_mplp_destroy(mplp);
	ave_depth = total_depth/(end - beg);
	return ave_depth;
}

void likelihood (bamFile fp,
				 bam_header_t* header,
				 bam_index_t* idx,
				double** transition, 
				 double** emission, 
				 char* ref, 
				 int32_t tid, 
				 int32_t window_beg,	// 0_based coordinate
				 int32_t region_beg,	// 0_based coordinate
				 int32_t region_end, 	// 0_based coordinate
				 int32_t size,
				 int32_t filter) {

	int32_t k, delet_count = 0;	// k is a relative coordinate within the window.
	for (k = region_beg - window_beg + 1; k < region_end - window_beg + 1; ++k) {	// change to 1_based coordinate
		if (delet_count > 0) {
			-- delet_count;
			continue;
		}
		if (ref[k - 1] == 'A' || ref[k - 1] == 'a' || ref[k - 1] == 'C' || ref[k - 1] == 'c' || ref[k - 1] == 'G' || 
		ref[k - 1] == 'g' || ref[k - 1] == 'T' || ref[k - 1] == 't') {

			int32_t beg = k + window_beg - 1 - size, end = k + window_beg - 1 + size;
			p_max* ref_allele = refp(emission, ref, k - 1);
			beg = beg < region_beg ? region_beg : beg;
			end = end > region_end ? region_end : end;
			
			/* Detect SNP. */
//			fprintf(stdout, "coordinate: %d\n", k + window_beg);
			if ((k + window_beg) == 2112653) {
				float test = base_read_depth(fp, idx, tid, k, beg, end);
				fprintf(stdout, "transition[%d][0] = %g, ref_allele->prob = %g, transition[%d][0] = %g, depth = %g\n", k - 1, transition[k - 1][0], ref_allele->prob, k, transition[k][0], test); 
			}

			if (transition[k - 1][0] >= 0.2 && ref_allele->prob <= 0.8 && transition[k][0] >= 0.2 && base_read_depth(fp, idx, tid, k, beg, end) > 5) {
				float qual = transition[k - 1][0] * transition[k][0];	// c*d
				double max;
				int8_t num;

				// Find out the base with highest emission probability.
				if (emission[k][1] > emission[k][2]) {
					max = emission[k][1];
					num = 1;
				} else {
					max = emission[k][2];
					num = 2;
				}
				if (max < emission[k][4]) {
					max = emission[k][4];
					num = 4;
				}
				if (max < emission[k][8]) {
					max = emission[k][8];
					num = 8;
				}
				if (max < emission[k][15]) {
					max = emission[k][15];
					num = 15;
				}

		//		if (max == 1) fprintf(stderr, "max: %g\tqual: %g\tnum: %d\n", max, qual, num);

				// Find out the 2nd max base.
				if (num != 15) {
					p_max* max2;
					switch (num) {
						case 1:
							max2 = bubble3(2, emission[k][2], 4, emission[k][4], 8, emission[k][8]);
							break;
						case 2:
							max2 = bubble3(1, emission[k][1], 4, emission[k][4], 8, emission[k][8]);
							break;
						case 4:
							max2 = bubble3(1, emission[k][1], 2, emission[k][2], 8, emission[k][8]);
							break;
						case 8:
							max2 = bubble3(1, emission[k][1], 2, emission[k][2], 4, emission[k][4]);
							break;
						default:
							fprintf(stderr, "The num is assigned wrongly.\n");
							exit (1);
					}
					
					if (num == ref_allele->num && max2->prob > 0.3) {	// max = ref allele
						char base = num2base(max2->num);
						qual = -4.343*log(1 - qual*max2->prob);
//						fprintf(stdout, "qual1: %g\n", qual);

						fprintf (stdout, "%s\t", header->target_name[tid]);
						fprintf (stdout, "%d\t.\t%c\t", k + window_beg, ref[k - 1]);
						fprintf (stdout, "%c\t%g\t", base, qual);
						if (filter == 0) fprintf (stdout, ".\t");
						else if (qual >= filter)	fprintf (stdout, "PASS\t");
						else fprintf (stdout, "q%d\t", filter);
						fprintf (stdout, "AF=%g\n", max2->prob);
					} else if (num != ref_allele->num){	// max != ref allele
						fprintf (stdout, "%s\t", header->target_name[tid]);
						fprintf (stdout, "%d\t.\t%c\t", k + window_beg, ref[k - 1]);
						qual = -4.343*log(1 - qual*max);
 //fprintf(stdout, "qual2: %g\n", qual);

						if (max2->prob > 0.3 && max2->num != ref_allele->num) {
							char base = num2base(num);
							char base2 = num2base(max2->num);
							fprintf(stdout, "%c,%c\t%g\t", base, base2, qual);
						}else{
							char base = num2base(num);
							fprintf(stdout, "%c\t%g\t", base, qual);
						} 
						if (filter == 0) fprintf (stdout, ".\t");
						else if (qual >= filter)	fprintf (stdout, "PASS\t");
						else fprintf (stdout, "q%d\t", filter);
						if (max2->prob > 0.3 && max2->num != ref_allele->num) fprintf (stdout, "AF=%g,AF=%g\n", max, max2->prob);
						else fprintf (stdout, "AF=%g\n", max);
					}
					free(max2);
				}
				free(ref_allele);
			}

			/* Detect insertion. */
			if (transition[k][1] > 0.3 && base_read_depth(fp, idx, tid, k, beg, end) > 5) {
				float qual = -4.343 * log(1 - transition[k][1]);
				float p = transition[k][1]/(transition[k][0] + transition[k][1]);
 
				fprintf (stdout, "%s\t%d\t.\t%c\t<I>\t%g\t", header->target_name[tid], k + window_beg, ref[k - 1], qual);
				if (filter == 0) fprintf (stdout, ".\t");
				else if (qual >= filter)	fprintf (stdout, "PASS\t");
				else fprintf (stdout, "q%d\t", filter);
				fprintf (stdout, "AF=%g\n", p);
			}

			/* Detect deletion. */	
			if (transition[k][2] > 0.3 && base_read_depth(fp, idx, tid, k, beg, end) > 5) {
				float diff = 0.3, qual;
				int32_t count1 = 1, count2 = 0, i;
				double path_p1 = transition[k][2], path_p2 = 0, path_ref = transition[k][0];

				// Record the 2 paths with highest probabilities.
				while (transition[k + count1][8] > transition[k + count1][7]) {
					float d = transition[k + count2][8] - transition[k + count2][7];
					if (d <= diff) {
						count2 = count1;
						path_p2 = path_p1*transition[k + count2][7];
						diff = d;
					}
					path_p1 *= transition[k + count2][8];
					++ count1;
				}				
				path_p1 *= transition[k + count1][7];
				fprintf (stdout, "%s\t%d\t.\t%c", header->target_name[tid], k + window_beg, ref[k - 1]);
				for (i = 0; i < count1; i ++) fprintf (stdout, "%c", ref[k + i]);
				fprintf (stdout, "\t%c", ref[k - 1]);

				for (i = 0; i < count2; ++i) {
					p_max* ref_allele = refp(emission, ref, k + i);
					path_ref *= (ref_allele->prob*transition[k + i][0]);
					free(ref_allele);
				}
				if (count2 > 0 && path_p2 > (path_ref*2)) {
					fprintf(stdout, ",%c", ref[k - 1]);
					for (i = k + count2; i < k + count1; i++) fprintf (stdout, "%c", ref[i]);
				} 
				qual = -4.343*log(1 - pow(path_p1, 1/count1));
				fprintf (stdout, "\t%g\t", qual);						
				if (filter == 0) fprintf (stdout, ".\t");
				else if (qual >= filter) fprintf (stdout, "PASS\t");
				else fprintf (stdout, "q%d\t", filter);
				if (count2 == 0 || (path_ref*2) >= path_p2) {
					float af = path_p1/(path_p1 + path_ref);
					fprintf (stdout, "AF=%g\n", af);
				} else {
					float total = path_p1 + path_p2;
					float af1 = path_p1/total;
					float af2 = path_p2/total;
					fprintf(stdout, "AF=%g,AF=%g\n", af1, af2);
				}
				delet_count = count1;
			}
		}
	}
}

