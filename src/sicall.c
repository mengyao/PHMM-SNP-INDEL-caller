/*
 * sicall.c: Calculate the probabilities of SNPs and INDELs. Make the calling.
 * Author: Mengyao Zhao
 * Create date: 2011-08-09
 * Contact: zhangmp@bc.edu
 * Last revise: 2014-06-11 
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
	for (iter = kh_begin(hi); iter != kh_end(hi); ++iter)
if (kh_exist(hi, iter)) fprintf(stderr, "pos: %d\tgenotype: %s\n", kh_key(hi, iter), kh_value(hi, iter).s);
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

void print_var (int32_t pos,
				int32_t filter,
				int32_t refa_begin,
				int32_t refa_end,
				char* ref_name,
				char* ref,
				char* var_allele1,
				char* var_allele2,
				double qual,
				double af1,
				double af2) {

	int32_t i;
	fprintf (stdout, "%s\t%d\t.\t", ref_name, pos);
	for (i = refa_begin; i < refa_end; ++i) fprintf(stdout, "%c", ref[i]);
	fprintf (stdout, "\t%s", var_allele1);
	if (af2 > 0) fprintf(stdout, ",%s", var_allele2);
	fprintf(stdout, "\t%g\t", qual);
	if (filter == 0) fprintf (stdout, ".\t");
	else if (qual >= filter)	fprintf (stdout, "PASS\t");
	else fprintf (stdout, "q%d\t", filter);
	fprintf (stdout, "AF=%g", af1);
	if (af2 > 0) fprintf(stdout, ",AF=%g", af2);
	fprintf(stdout, "\n");
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

	int32_t k, jump_count = 0, k_beg = region_beg > window_beg ? region_beg - window_beg + 1 : 1;	// k is a relative coordinate within the window.
	for (k = k_beg; k < region_end - window_beg + 1; ++k) {	// change to 1_based coordinate
		if (jump_count > 0) {
			-- jump_count;
			continue;
		}
		if (ref[k - 1] == 'A' || ref[k - 1] == 'a' || ref[k - 1] == 'C' || ref[k - 1] == 'c' || ref[k - 1] == 'G' || ref[k - 1] == 'g' || ref[k - 1] == 'T' || ref[k - 1] == 't') {
	//fprintf(stderr, "ref[%d]: %c\n", k - 1, ref[k - 1]);

			int32_t beg = k - 1 - size, end = k - 1 + size;
			char* var_allele1 = malloc(100*sizeof(char));
			char* var_allele2 = malloc(100*sizeof(char));
	
			p_max* ref_allele = refp(emission, ref, k - 1);
			beg = beg < 0 ? 0 : beg;
			end = end > region_end - window_beg ? region_end - window_beg : end;
			var_allele2[0] = '\0';

			/* Detect SNP. */
			if (transition[k - 1][0] >= 0.2 && ref_allele->prob <= 0.8 && transition[k][0] >= 0.2) {
				p_cov c = cov(cinfo, beg, end);
				double qual = transition[k - 1][0] * transition[k][0];	// c*d
				p_max* max = max1and2(emission[k][1], emission[k][2], emission[k][4], emission[k][8], emission[k][15]);
				if (max[0].base != 'N' && c.ave_depth > 5 && c.map_qual >= 10) {
					if (max[0].base == ref_allele->base && max[1].prob > 0.3) {	// max = ref allele
						var_allele1[0] = max[1].base;
						var_allele1[1] = '\0';
						qual = -4.343*log(1 - qual*max[1].prob);
						print_var (k + window_beg, filter, k - 1, k, header->target_name[tid], ref, var_allele1, var_allele2, qual, max[1].prob, 0);
					} else if (max[0].base != ref_allele->base){	// max != ref allele, this is where the error snps come from 
						double af2 = 0;
						var_allele1[0] = max[0].base;
						var_allele1[1] = '\0';
						qual = -4.343*log(1 - qual*max[0].prob);
						if (max[1].prob > 0.3 && max[1].base != ref_allele->base) {
							var_allele2[0] = max[1].base;
							var_allele2[1] = '\0';
							af2 = max[1].prob;
						}
						print_var (k + window_beg, filter, k - 1, k, header->target_name[tid], ref, var_allele1, var_allele2, qual, max[0].prob, af2);
					}
				}//
				free (max);
				free(ref_allele);
			}

			/* Detect insertion. */
			if (transition[k][1] > 0.3 && transition[k][2] < 0.3 && k < region_end - window_beg) {
				p_cov c = cov(cinfo, beg, end);
				p_haplotype* haplo = haplotype_construct(hi, hm, hd, 1, k);
				if (haplo && haplo->haplotype1[0] != 'N' && c.ave_depth > 5 && c.map_qual >= 10) {
					double qual, p;
					int32_t i = k + 1, indel_dis = 0;
					while (ref[i] == haplo->haplotype1[0]) {
						++ indel_dis;
						++ i;
					}
					if (strlen(haplo->haplotype1) == 1 && transition[i][2] > 0.3) {	// SNP presented as INDEL
						int32_t j = i + 1, delet_len;
						double af1 = haplo->count1/c.ave_depth;
						p_haplotype* haplod = haplotype_construct(hi, hm, hd, 2, i);
						while (ref[j] == ref[i]) {
							++ indel_dis;
							++ j;
						}
						if (haplod) {
							delet_len = strlen(haplod->haplotype1);
							haplotype_destroy (haplod);
						} else delet_len = 1;
						if (delet_len > j - i) indel_dis += delet_len - j + i;

						p = (i - k)*transition[k][1];
						p = p > 1 ? 1 : p;
						qual = -4.343 * log(1 - p);
						af1 = af1 > 1 ? 1 : af1;
						print_var (i + window_beg + 1, filter, i, i + delet_len, header->target_name[tid], ref, haplo->haplotype1, var_allele2, qual, af1, 0);
						jump_count = indel_dis;
					} else { 
fprintf(stderr, "t[%d]: %g\tcount: %d\n", k, transition[k][1], haplo->count1);
						double af1, af2;
						char refa[] = {ref[k - 1], '\0'};
						p = haplo->count1/c.ave_depth;
						p = p > 1 ? 1 : p;
						//qual = -4.343 * log(1 - transition[k][1]);
						qual = -4.343 * log(1 - p);
						strcpy(var_allele1, refa);
						if (haplo->count2 == 0) {
							af1 = p;
							af2 = 0;
						}else {
							af1 = p*haplo->count1/(haplo->count1 + haplo->count2);
							if (haplo->count2 > 5) {
								strcpy(var_allele2, refa);
								af2 = p*haplo->count2/(haplo->count1 + haplo->count2);
								strcat(var_allele2, haplo->haplotype2);
							} else af2 = 0;
						}
						print_var (k + window_beg, filter, k - 1, k, header->target_name[tid], ref, strcat(var_allele1, haplo->haplotype1), var_allele2, qual, af1, af2);
					}
					haplotype_destroy(haplo);
				}//
			}

			/* Detect deletion. */
			// homopolymer deletion
			if (k + 1 <= strlen(ref) && ref[k + 1] == ref[k]) {	// ref: 0-based
				p_cov c = cov(cinfo, beg, end);	// cov return read depth and mapping quality
				if (c.ave_depth > 5 && c.map_qual >= 10) {
					int32_t mer_len = 1, delet_len = 0, i, l = 0, pos = 0, seg_count = 0, skip_len = 0;
					double t = 0, p = 1, af, afs = 0;
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
									pos = k + i;
									skip_len = l + i;
									break;
								} 
							} else haplotype_destroy (haplo);
						}
					}
					if (haplo && l > 0 && pos > 0 && haplo->count1 > 4) {	// deletion containing bases after the homopolymer
						double qual, pl, af1 = afs/seg_count;
						pl = transition[pos][2];
						for (i = 1; i < l; ++i) pl *= transition[pos + i][8];
						pl *= transition[pos + i + 1][7];
						pl = pow(pl, 1/l);
						p = p > pl ? p : pl;
						p = p > 1 ? 1 : p;
						qual = -4.343 * log(1 - p);
						var_allele1[0] = ref[pos - 1];
						var_allele1[1] = '\0';
						print_var (pos + window_beg, filter, pos - 1, pos + l, header->target_name[tid], ref, var_allele1, var_allele2, qual, af1, 0);
						jump_count = skip_len >= mer_len ? skip_len : (mer_len - 1);
					} else if (t > 0.3 && delet_len > 0 && delet_len <= mer_len && transition[k][1] < 0.3) {

						int32_t indel_dis = mer_len, mnp = 0, end;
						double qual = -4.343 * log(1 - p), af1 = afs/seg_count;
						p_haplotype* haploi = 0;					
						end = k + mer_len;
						while (ref[end] == ref[k + mer_len]) ++end;
						for (i = k + 1; i <= end; ++i) {
							if (transition[i][1] > 0.2) {
								haploi = haplotype_construct(hi, hm, hd, 1, i);
								if (haploi && haploi->count1 > 4 && strlen(haploi->haplotype1) == delet_len) {
									mnp = 1;
									i = end + 1;
								}
							}
							++indel_dis;
						}	
						if (transition[i][1] > 0.2)	haploi = haplotype_construct(hi, hm, hd, 1, i);
						if (haploi && haploi->count1 > 4 && strlen(haploi->haplotype1) == delet_len) mnp = 1;
						if (mnp) {	// MNP is called as INDEL
							print_var (k + mer_len + window_beg, filter, k, k + delet_len, header->target_name[tid], ref, haploi->haplotype1, var_allele2, qual, af1, 0);
							jump_count = indel_dis;
						} else {
							var_allele1[0] = ref[k - 1];
							var_allele1[1] = '\0';
							print_var (k + window_beg, filter, k - 1, k + delet_len, header->target_name[tid], ref, var_allele1, var_allele2, qual, af1, 0);
							jump_count = delet_len >= mer_len ? delet_len : (mer_len - 1);
						}
						if (haploi) haplotype_destroy(haploi);
					}
					if (haplo && pos > 0) haplotype_destroy(haplo);
				}//
			} else if (transition[k][2] > 0.3) {	// transition: 1-based
				p_cov c = cov(cinfo, beg, end);
//fprintf(stderr, "k: %d************************************************\nc.ave_depth: %g\tc.map_qual: %g\n", k, c.ave_depth, c.map_qual);
				if (c.ave_depth > 4.5 && c.map_qual >= 10) {
					double diff = 0.3, qual, total, af1, af2, p;
					int32_t count1, count2 = 0, i;
					double path_p1 = transition[k][2], path_p2 = 0;//, path_ref = transition[k][0];

					if (transition[k + 1][2] > 0.3) {	// symitry deletion region
//fprintf(stderr, "here\n");
						p_haplotype* haplo; 
						for (i = 1; i <= 2; ++i) {
							haplo = haplotype_construct(hi, hm, hd, 2, k + i);
							if (haplo && haplo->count1 > 4) {
								af1 = haplo->count1/c.ave_depth;
								if (af1 > 0.3) {
									af1 = af1 > 1 ? 1 : af1;
									qual = -4.343*log(1 - af1);
									var_allele1[0] = ref[k + i - 2];
									var_allele1[1] = '\0';
									count1 = strlen(haplo->haplotype1);
									print_var (k + i + window_beg - 1, filter, k + i - 2, k + count1 + i - 1, header->target_name[tid], ref, var_allele1, var_allele2, qual, af1, 0);
									jump_count = i + count1 - 1;
									i = 2;
								}
							}
							if (haplo) haplotype_destroy(haplo);
						}
					} else {
						p_haplotype* haplo = haplotype_construct(hi, hm, hd, 2, k + 1);
//if (haplo) fprintf(stderr, "haplo->count: %d\n", haplo->count1);
						if (haplo && haplo->count1/c.ave_depth > 0.3) {
						// Record the 2 paths with highest probabilities.
							count1 = 1;
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
							var_allele1[0] = ref[k - 1];
							var_allele1[1] = '\0';
							p = haplo->count1/c.ave_depth;
							p = p > 1 ? 1 : p;
							qual = -4.343*log(1 - p);

							total = path_p1 + path_p2;
							af2 = path_p2/total;

							if (af2 > 0.3 && count2 > 0) { 
								int32_t n = 1;
								af1 = path_p1/total;
								var_allele2[0] = ref[k - 1];
								for (i = k + count2; i < k + count1; i++) var_allele2[n ++] = ref[i];
								var_allele2[n] = '\0';	
								print_var (k + window_beg, filter, k - 1, k + count1, header->target_name[tid], ref, var_allele1, var_allele2, qual, af1, af2);
							} else print_var (k + window_beg, filter, k - 1, k + count1, header->target_name[tid], ref, var_allele1, var_allele2, qual, p, 0);
							jump_count = count1;
						}
						if (haplo) haplotype_destroy(haplo);
					}
				}//
			}
			free (var_allele2);
			free (var_allele1);
		}
	}
}

