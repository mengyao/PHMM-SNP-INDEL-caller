/*
 * sicall.c: Calculate the probabilities of SNPs and INDELs. Make the calling.
 * Author: Mengyao Zhao
 * Create date: 2011-08-09
 * Contact: zhangmp@bc.edu
 * Last revise: 2014-04-11 
 */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "khash.h"
#include "sicall.h"
#include "viterbi.h"

#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

#ifndef KHASH
#define KHASH
KHASH_MAP_INIT_INT(insert, kstring_t)
KHASH_MAP_INIT_INT(mnp, kstring_t)
KHASH_MAP_INIT_INT(delet, kstring_t)
#endif
KHASH_MAP_INIT_STR(count, int32_t)

typedef struct {
	char base;
	double prob;
} p_max;

typedef struct {
	char* haplotype1;
	char* haplotype2;
	int32_t count1;
	int32_t count2;
} p_haplotype;

typedef struct {
	float ave_depth;
	float map_qual;
} p_cov;

// Return the base and emission probability of ref_allele.
p_max* refp (double** emission, char* ref, int32_t k) {
	p_max* ref_allele = (p_max*)malloc(sizeof(p_max));
	switch (ref[k]) {
		case 'A':
		case 'a':
			ref_allele->prob = emission[k + 1][1];
			ref_allele->base = 'A';
			break;
		case 'C':
		case 'c':
			ref_allele->prob = emission[k + 1][2];
			ref_allele->base = 'C';
			break;
		case 'G':
		case 'g':
			ref_allele->prob = emission[k + 1][4];
			ref_allele->base = 'G';
			break;
		case 'T':
		case 't':
			ref_allele->prob = emission[k + 1][8];
			ref_allele->base = 'T';
			break;
		default:
			fprintf(stderr, "Wrong reference sequence. \n");
			exit (1);
			break;
	}
	return ref_allele;
}

p_max bubble3 (char c1, double p1, char c2, double p2, char c3, double p3) {
	p_max m;// = (p_max*)malloc(sizeof(p_max));
	m.base = c1;
	m.prob = p1;
	if (p2 > m.prob) {
		m.base = c2;
		m.prob = p2;
	}
	if (p3 > m.prob) {
		m.base = c3;
		m.prob = p3;
	}
	return m;
}

p_max* max1and2 (double p1, double p2, double p3, double p4, double p5) {
	p_max *m = (p_max*)malloc(2*sizeof(p_max));
	// Find out the base with highest emission probability.
	if (p1 > p2) {
		m[0].prob = p1;
		m[0].base = 'A';
	} else {
		m[0].prob = p2;
		m[0].base = 'C';
	}
	if (m[0].prob < p3) {
		m[0].prob = p3;
		m[0].base = 'G';
	}
	if (m[0].prob < p4) {
		m[0].prob = p4;
		m[0].base = 'T';
	}
	if (m[0].prob < p5) {
		m[0].prob = p5;
		m[0].base = 'N';
	}

	// Find out the 2nd max base.
	if (m[0].base != 'N') {
		switch (m[0].base) {
			case 'A':
				m[1] = bubble3('C', p2, 'G', p3, 'T', p4);
				break;
			case 'C':
				m[1] = bubble3('A', p1, 'G', p3, 'T', p4);
				break;
			case 'G':
				m[1] = bubble3('A', p1, 'C', p2, 'T', p4);
				break;
			case 'T':
				m[1] = bubble3('A', p1, 'C', p2, 'G', p3);
				break;
			default:
				fprintf(stderr, "The max base is wrong.\n");
				exit (1);
		}
	}
	return m;
}

//filter based on read depth and mapping quality of the region
p_cov cov(p_info* cinfo, int32_t beg, int32_t end) {
	float sum = 0, qual = 0;
	int32_t i;
	p_cov r;
	for (i = beg; i <= end; ++i) {
		sum += cinfo[i].depth;
		if (cinfo[i].depth > 0) qual += (cinfo[i].mqual_sum/cinfo[i].depth);
		else qual += 0;
	}
	sum /= (end - beg + 1);
	qual /=(end - beg + 1);
	r.ave_depth = sum;
	r.map_qual = qual;
	return r;
}

// Combine repeat genotypes and prepare for printing.
p_haplotype* haplotype_construct (khash_t(insert) *hi,
								  khash_t(mnp) *hm,
								  khash_t(delet) *hd,
							   	  int32_t type,	// 0: mnp, 1: insert, 2: delet
								  int32_t pos){

	khiter_t iter;
	int32_t i, len = 128;
	char* genotype;
	int32_t c = 0, total_len;
	int ret;
	khash_t(count) *hc = kh_init(count);
	khiter_t ic;
/*
fprintf(stderr, "haplotype construct\n");
	for (iter = kh_begin(hd); iter != kh_end(hd); ++iter)
if (kh_exist(hd, iter)) fprintf(stderr, "pos: %d\tgenotype: %s\n", kh_key(hd, iter), kh_value(hd, iter).s);
fprintf(stderr, "type: %d\tpos: %d\n", type, pos);
*/

	if (type == 0) {	
		iter = kh_get(mnp, hm, pos);	//
		if (iter == kh_end(hm)) return 0;	//
	} else if (type == 1) {
		iter = kh_get(insert, hi, pos);	//
		if (iter == kh_end(hi)) return 0;	//
	}else if (type == 2) {
		iter = kh_get(delet, hd, pos);	//
		if (iter == kh_end(hd)) return 0;	//
	}

	char* key = (char*)malloc(len*sizeof(char));
	p_haplotype* h = (p_haplotype*)malloc(sizeof(p_haplotype));
	
	if (type == 0) genotype = kh_value(hm, iter).s;	//
	else if (type == 1) genotype = kh_value(hi, iter).s;	//
	else if (type == 2) genotype = kh_value(hd, iter).s;	//
	else {
		fprintf(stdout, "The parameter type of the function haplotype_construct can only be assigned as 0, 1, or 2.\n");
		exit(1);
	}

	total_len = strlen(genotype);
	for (i = 0; i < total_len; ++i) {
		if (genotype[i] == ',' || i == (total_len - 1)) {
			key[c] = '\0';
			c = 0;
			ic = kh_put(count, hc, key, &ret);
			if (ret == 0) kh_value(hc, ic) = kh_value(hc, ic) + 1;	// The key exist.
			else kh_value(hc, ic) = 1;	// The key doesn't exist.
			key = (char*)malloc(len*sizeof(char));
		} else {
			if (c + 2 >= len) {
				++len;
				kroundup32(len);
				key = realloc(key, len * sizeof(int32_t));	
			}
			key[c++] = genotype[i]; 
		}
	}

	h->count1 = 0;
	h->haplotype1 = (char*)malloc(len*sizeof(char));
	for(ic = kh_begin(hc); ic != kh_end(hc); ++ic) {
		if (kh_exist(hc, ic) && kh_value(hc, ic) > h->count1) {
			h->count1 = kh_value(hc, ic);
			strcpy(h->haplotype1, kh_key(hc, ic));
		}
	}
	ic = kh_get(count, hc, h->haplotype1);
	kh_del(count, hc, ic);
	h->count2 = 0;
	h->haplotype2 = (char*)malloc(len*sizeof(char));
	for(ic = kh_begin(hc); ic != kh_end(hc); ++ic) {
		if (kh_exist(hc, ic) && kh_value(hc, ic) > h->count2) {
			h->count2 = kh_value(hc, ic);
			strcpy(h->haplotype2, kh_key(hc, ic));
		}
	}

	for (ic = kh_begin(hc); ic != kh_end(hc); ++ic) {
		if (!kh_exist(hc,ic)) continue;
		key = (char*)kh_key(hc, ic);
		free(key);
	}
	kh_destroy(count, hc);

	return h;
}

void haplotype_destroy (p_haplotype* hapo) {
	free(hapo->haplotype1);
	free(hapo->haplotype2);
	free(hapo);
}

void likelihood (bam_header_t* header,
				double** transition, 
				 double** emission, 
				 char* ref,
				p_info* cinfo, 
				 int32_t tid, 
				 int32_t window_beg,	// 0_based coordinate
				 int32_t region_beg,	// 0_based coordinate
				 int32_t region_end, 	// 0_based coordinate
				 int32_t size,
				 int32_t filter,
				khash_t(insert) *hi,
				khash_t(mnp) *hm,
				khash_t(delet) *hd) {

	int32_t k, delet_count = 0, k_beg = region_beg > window_beg ? region_beg - window_beg + 1 : 1;	// k is a relative coordinate within the window.
//fprintf(stderr, "region_beg: %d\twindow_beg: %d\n", region_beg, window_beg);
	for (k = k_beg; k < region_end - window_beg + 1; ++k) {	// change to 1_based coordinate
		if (delet_count > 0) {
			-- delet_count;
			continue;
		}
		if (ref[k - 1] == 'A' || ref[k - 1] == 'a' || ref[k - 1] == 'C' || ref[k - 1] == 'c' || ref[k - 1] == 'G' || ref[k - 1] == 'g' || ref[k - 1] == 'T' || ref[k - 1] == 't') {
//	fprintf(stderr, "ref[%d]: %c\n", k - 1, ref[k - 1]);

			int32_t beg = k - 1 - size, end = k - 1 + size;
			p_max* ref_allele = refp(emission, ref, k - 1);
			beg = beg < 0 ? 0 : beg;
			end = end > region_end - window_beg ? region_end - window_beg : end;

			/* Detect SNP. */
			if (transition[k - 1][0] >= 0.2 && ref_allele->prob <= 0.8 && transition[k][0] >= 0.2) {
				p_cov c = cov(cinfo, beg, end);
				float qual = transition[k - 1][0] * transition[k][0];	// c*d
				p_max* max = max1and2(emission[k][1], emission[k][2], emission[k][4], emission[k][8], emission[k][15]);
				if (max[0].base != 'N' && c.ave_depth > 5 && c.map_qual >= 10) {
					if (max[0].base == ref_allele->base && max[1].prob > 0.3) {	// max = ref allele
						qual = -4.343*log(1 - qual*max[1].prob);
						fprintf (stdout, "%s\t", header->target_name[tid]);
						fprintf (stdout, "%d\t.\t%c\t", k + window_beg, ref[k - 1]);
						fprintf (stdout, "%c\t%g\t", max[1].base, qual);
						if (filter == 0) fprintf (stdout, ".\t");
						else if (qual >= filter)	fprintf (stdout, "PASS\t");
						else fprintf (stdout, "q%d\t", filter);
						fprintf (stdout, "AF=%g\n", max[1].prob);
					} else if (max[0].base != ref_allele->base){	// max != ref allele, this is where the error snps come from 
						fprintf (stdout, "%s\t", header->target_name[tid]);
						fprintf (stdout, "%d\t.\t%c\t", k + window_beg, ref[k - 1]);
						qual = -4.343*log(1 - qual*max[0].prob);
						if (max[1].prob > 0.3 && max[1].base != ref_allele->base) fprintf(stdout, "%c,%c\t%g\t", max[0].base, max[1].base, qual);
						else fprintf(stdout, "%c\t%g\t", max[0].base, qual);
						if (filter == 0) fprintf (stdout, ".\t");
						else if (qual >= filter)	fprintf (stdout, "PASS\t");
						else fprintf (stdout, "q%d\t", filter);
						if (max[1].prob > 0.3 && max[1].base != ref_allele->base) fprintf (stdout, "AF=%g,AF=%g\n", max[0].prob, max[1].prob);
						else fprintf (stdout, "AF=%g\n", max[0].prob);
					}
				}
				free (max);
				free(ref_allele);
			}

			/* Detect insertion. */
			if (transition[k][1] > 0.3 && transition[k][2] < 0.3) {
				p_cov c = cov(cinfo, beg, end);
				p_haplotype* haplo = haplotype_construct(hi, hm, hd, 1, k);
				if (haplo && c.ave_depth > 5 && c.map_qual >= 10) {
					float qual = -4.343 * log(1 - transition[k][1]);
					float p = transition[k][1]/(transition[k][0] + transition[k][1]);
					if (strlen(haplo->haplotype1) == 1) {
						p_max* max = max1and2(emission[k][0], emission[k][3], emission[k][5], emission[k][9], emission[k][14]);
						fprintf (stdout, "%s\t%d\t.\t%c\t%c%c", header->target_name[tid], k + window_beg, ref[k - 1], ref[k - 1], max[0].base);
						if(max[1].prob > 0.3) fprintf(stdout, ",%c%c", ref[k - 1], max[1].base);
						fprintf(stdout, "\t%g\t", qual);
						if (filter == 0) fprintf (stdout, ".\t");
						else if (qual >= filter)	fprintf (stdout, "PASS\t");
						else fprintf (stdout, "q%d\t", filter);
						if (max[1].prob == 0)fprintf (stdout, "AF=%g\n", p);
						else {
							float p1 = p*max[0].prob/(max[0].prob + max[1].prob);
							fprintf(stdout, "AF=%g", p1);
							if (max[1].prob > 0.3) {
								float p2 = p*max[1].prob/(max[0].prob + max[1].prob);
								fprintf(stdout, ",%g", p2);
							}
							fprintf(stdout, "\n");
						}
					} else {
						fprintf (stdout, "%s\t%d\t.\t%c\t%c%s", header->target_name[tid], k + window_beg, ref[k - 1], ref[k - 1], haplo->haplotype1);
						if(haplo->count2 > 5) fprintf(stdout, ",%c%s", ref[k - 1], haplo->haplotype2);
						fprintf(stdout, "\t%g\t", qual);
						if (filter == 0) fprintf (stdout, ".\t");
						else if (qual >= filter)	fprintf (stdout, "PASS\t");
						else fprintf (stdout, "q%d\t", filter);
						if (haplo->count2 == 0)fprintf (stdout, "AF=%g\n", p);
						else {
							float p1 = p*haplo->count1/(haplo->count1 + haplo->count2);
							fprintf(stdout, "AF=%g", p1);
							if (haplo->count2 > 5) {
								float p2 = p*haplo->count2/(haplo->count1 + haplo->count2);
								fprintf(stdout, ",%g", p2);
							}
							fprintf(stdout, "\n");
						}
					}
					haplotype_destroy(haplo);
				}//
			}

			/* Detect deletion. */
			// homopolymer deletion
			if (k + 2 <= strlen(ref) && ref[k + 1] == ref[k] && ref[k + 2] == ref[k]) {	// ref: 0-based
				p_cov c = cov(cinfo, beg, end);	// cov return read depth and mapping quality
				if (c.ave_depth > 5 && c.map_qual >= 10) {
					int32_t mer_len = 1, delet_len = 0, i, l = 0, pos = 0, seg_count = 0, skip_len = 0;
					float t = 0, p = 1, af, afs = 0;
					p_haplotype* haplo;
					while (ref[k + mer_len] == ref[k]) ++ mer_len;
					for (i = 0; i < mer_len; ++i) {
						haplo = haplotype_construct(hi, hm, hd, 2, k + i + 1);
						if (haplo) {
							af = haplo->count1/c.ave_depth;
							af = af > 1 ? 1 : af;
							if (af > 0.3) {
								int32_t j; 
								l = (int32_t)strlen(haplo->haplotype1);
								afs += af;
								++seg_count;
								if (l + i <= mer_len) {
									t += transition[k + i][2];
									delet_len += l;
									haplotype_destroy (haplo);
									p *= transition[k + i][2];
									for (j = 1; j < l; ++j) p *= transition[k + i + j][8];
									p *= transition[k + i + j][7];
									p = pow(p, 1/(l + 1));
								} else {	// deletion containing bases after the homopolymer
									pos = k + i + 1;
									skip_len = l + i;
									break;
								} 
							} else haplotype_destroy (haplo);
						}
					}
					if (haplo && l > 0 && pos > 0) {	// deletion containing bases after the homopolymer
						float qual, pl;
						pl = transition[pos][2];
						for (i = 1; i < l; ++i) pl *= transition[pos + i][8];
						pl *= transition[pos + i + 1][7];
						pl = pow(pl, 1/l);
						p = p > pl ? p : pl;
						qual = -4.343 * log(1 - p);
						fprintf (stdout, "%s\t%d\t.\t%c", header->target_name[tid], pos - delet_len + window_beg, ref[pos - delet_len - 1]);
						for (i = 0; i < l + delet_len; ++i) fprintf(stdout, "%c", ref[pos - delet_len + i]);
						fprintf(stdout, "\t%c\t%g\t", ref[pos - delet_len - 1], qual);
						if (filter == 0) fprintf (stdout, ".\t");
						else if (qual >= filter)	fprintf (stdout, "PASS\t");
						else fprintf (stdout, "q%d\t", filter);
						fprintf(stdout, "AF=%g\n", afs/seg_count);
						haplotype_destroy (haplo);
					} else if (t > 0.3 && delet_len > 0 && delet_len <= mer_len && transition[k][1] < 0.3) {
						float qual = -4.343 * log(1 - p);
						fprintf (stdout, "%s\t%d\t.\t%c", header->target_name[tid], k + window_beg, ref[k - 1]);
						for (i = 0; i < delet_len; ++i) fprintf(stdout, "%c", ref[k + i]);
						fprintf(stdout, "\t%c\t%g\t", ref[k - 1], qual);
						if (filter == 0) fprintf (stdout, ".\t");
						else if (qual >= filter)	fprintf (stdout, "PASS\t");
						else fprintf (stdout, "q%d\t", filter);
						fprintf(stdout, "AF=%g\n", afs/seg_count);
						skip_len = delet_len;
					}
					skip_len = skip_len >= mer_len ? skip_len : (mer_len - 1);
					delet_count = skip_len;
				}//
			} else if (transition[k][2] > 0.3) {	// transition: 1-based
				p_cov c = cov(cinfo, beg, end);
				if (c.ave_depth > 5 && c.map_qual >= 10) {
					float diff = 0.3, qual, total, af1, af2;
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
					fprintf (stdout, "\t%c", ref[k - 1]);	// the reference base before deletion

					for (i = 0; i < count2; ++i) {
						p_max* ref_allele = refp(emission, ref, k + i);
						path_ref *= (ref_allele->prob*transition[k + i][0]);
						free(ref_allele);
					}

					total = path_p1 + path_p2;
					af1 = path_p1/total;
					af2 = path_p2/total;
					if (af2 > 0.01 && count2 > 0 && path_p2 > (path_ref*2)) {
						fprintf(stdout, ",%c", ref[k - 1]);
						for (i = k + count2; i < k + count1; i++) fprintf (stdout, "%c", ref[i]);
					} 
					qual = -4.343*log(1 - pow(path_p1, 1/count1));
					fprintf (stdout, "\t%g\t", qual);						
					if (filter == 0) fprintf (stdout, ".\t");
					else if (qual >= filter) fprintf (stdout, "PASS\t");
					else fprintf (stdout, "q%d\t", filter);
					
					if (af2 > 0.01 && count2 > 0 && path_p2 > (path_ref*2)) fprintf(stdout, "AF=%g,AF=%g\n", af1, af2);
					else {
						float af = path_p1/(path_p1 + path_ref);
						fprintf (stdout, "AF=%g\n", af);
					}
					delet_count = count1;
				}//
			}
		}
	}
}

