/*
 * sicall.c: Calculate the probabilities of SNPs and INDELs. Make the calling.
 * Author: Mengyao Zhao
 * Create date: 2011-08-09
 * Contact: zhangmp@bc.edu
 * Last revise: 2013-01-29 
 */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "sicall.h"

typedef struct {
	int8_t num;
	double prob;
} p_max;

p_max* refp (double** emission, char* ref, int32_t k) {
	p_max* ref;
	switch (ref[k]) {
		case 'A':
		case 'a':
			ref->prob = emission[k + 1][1];
			ref->num = 1;
			break;
		case 'C':
		case 'c':
			ref->prob = emission[k + 1][2];
			ref->num = 2;
			break;
		case 'G':
		case 'g':
			ref->prob = emission[k + 1][4];
			ref->num = 4;
			break;
		case 'T':
		case 't':
			ref->prob = emission[k + 1][8];
			ref->num = 8;
			break;
		default:
			fprintf(stderr, "Wrong reference sequence. \n");
			exit (1);
			break;
	}
	return ref;
}

p_max* bubble3 (int8_t n1, double p1, int8_t n2, double p2, int8_t n3, double p3) {
//	max *m = (max*)malloc(sizeof(max));
	p_max *m;
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

void likelihood (double** transition, 
				 double** emission, 
				 char* ref, 
				 char* ref_name, 
				 int32_t window_beg,	// 0_based coordinate
				 int32_t region_beg,	// 0_based coordinate
				 int32_t region_end, 	// 0_based coordinate
				 int32_t filter) {

	int32_t k;
	for (k = region_beg - window_beg + 1; k < region_end - window_beg + 1; ++k) {	// change to 1_based coordinate
		if (ref[k - 1] == 'A' || ref[k - 1] == 'a' || ref[k - 1] == 'C' || ref[k - 1] == 'c' || ref[k - 1] == 'G' || 
		ref[k - 1] == 'g' || ref[k - 1] == 'T' || ref[k - 1] == 't') {

			/* Detect SNP. */
			p_max* ref = refp(emission, ref, k - 1);
		//	double e;
		//	int8_t ref_num;
			
			if (transition[k - 1][0] >= 0.2 && ref->prob <= 0.8 && transition[k][0] >= 0.2) {
				float qual = transition[k - 1][0] * transition[k][0];	// c*d
				double max;
				int8_t num;
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
					
					if (num == ref->num && max2->prob > 0.3) {	// max = ref allele
						char base = num2base(max2->num);
						qual = -4.343*log(1 - qual*max2->prob);
						fprintf (stdout, "%s\t", ref_name);
						fprintf (stdout, "%d\t.\t%c\t", k + window_beg, ref[k - 1]);
						fprintf (stdout, "%c\t%f\t", base, qual);
						if (filter == 0) fprintf (stdout, ".\t");
						else if (qual >= filter)	fprintf (stdout, "PASS\t");
						else fprintf (stdout, "q%d\t", filter);
						fprintf (stdout, "AF=%f\n", max2->prob);
					} else {	// max != ref allele
						fprintf (stdout, "%s\t", ref_name);
						fprintf (stdout, "%d\t.\t%c\t", k + window_beg, ref[k - 1]);
						if (max2->prob > 0.3 && max2->num != ref->num) {
							char base = num2base(num);
							char base2 = num2base(max2->num);
							fprintf(stdout, "%c,%c\t", base, base2);
						}else{
							char base = num2base(num);
							fprintf(stdout, "%c\t", base);
						} 
						qual = -4.343*log(1 - qual*max);
						if (filter == 0) fprintf (stdout, ".\t");
						else if (qual >= filter)	fprintf (stdout, "PASS\t");
						else fprintf (stdout, "q%d\t", filter);
						if (max2->prob > 0.3 && max2->num != ref->num) fprintf (stdout, "AF=%f,AF=%f\n", max, max2->prob);
						else fprintf (stdout, "AF=%f\n", max);
					}
				}

			/* Detect insertion. */
			if (transition[k][1] > 0.3) {
				float qual = -4.343 * log(1 - transition[k][1]);
				float p = transition[k][1]/transition[k][0]; 
				fprintf (stdout, "%s\t%d\t.\t%c\t<I>\t%f\t", ref_name, k + window_beg, ref[k - 1], qual);
				if (filter == 0) fprintf (stdout, ".\t");
				else if (qual >= filter)	fprintf (stdout, "PASS\t");
				else fprintf (stdout, "q%d\t", filter);
				fprintf (stdout, "AF=%f\n", p);
			}

			/* Detect deletion. */	
			if (transition[k][2] > 0.3) {
				// Record the 2 paths with highest probabilities.
				double path_p = transition[k][2];
				double path_p2 = transition[k][2];
				int32_t count = 1, count2 = 1, flag = 0;//FIXME
				while (transition[k + count2][7] < transition[k + count2][8]) {
					if (flag == 0) {
						path_p *= transition[k + count2][8];
						path_p2 *= transition[k + count2][7];
					} else path_p2 *= transition[k + count2][8];
					if (path_p < 0.1) path_p = 0;
					else if (path_p2 < 0.1 && flag == 0) path_p2 = path_p;
					else if (path_p2 < 0.1 && flag == 1) path_p2 = 0;
					else flag = 1;	/* path_p is determined */
					if (flag == 0) count ++;
					count2 ++;
				}
				if (flag == 0) {
					path_p *= transition[k + count2][7];
					path_p2 *= transition[k + count2][8];
				} else path_p2 *= transition[k + count2][7];
				if (path_p < 0.1) path_p = 0;
				else if (path_p2 < 0.1) path_p2 = 0;
				else if (flag == 0) {
					count2 ++;
					while (transition[k + count2][7] < transition[k + count2][8]) {
						path_p2 *= transition[k + count2][8];
						if (path_p2 < 0.1) {
							path_p2 = 0;
							break;
						}
						count2 ++;
					}
					path_p2 *= transition[k + count2][7];
					if (path_p2 < 0.1) path_p2 = 0;
				}

				if (path_p >= 0.1) {
					int32_t count_max = path_p2 >= 0.1 ? count2 : count;
					int32_t i; 
					double t = transition[k][0];	// probability of the main path
					int32_t count_p = path_p2 > path_p ? count2 : count;
					double p = path_p2 > path_p ? path_p2 : path_p;
					for (i = 1; i <= count_p; ++i) t *= (emission[k + i][]*transition[k + i][0]);
				//	float qual = -4.343 * log(t/(t + p));
					float qual = -4.343*log(1 - p/count_p);
					fprintf (stdout, "%s\t%d\t.\t%c", ref_name, k + window_beg, ref[k - 1]);
					for (i = 0; i < count_max; i ++) fprintf (stdout, "%c", ref[k + i]);
					fprintf (stdout, "\t%c", ref[k - 1]);
					for (i = k + count; i < k + count2; i++) fprintf (stdout, "%c", ref[i]);

					if (path_p2 >= 0.1 && count2 > count) fprintf (stdout, ",%c",ref[k - 1]);
					fprintf (stdout, "\t%f\t", qual);						
					if (filter == 0) fprintf (stdout, ".\t");
					else if (qual >= filter)	fprintf (stdout, "PASS\t");
					else fprintf (stdout, "q%d\t", filter);
					fprintf (stdout, "AF=%f", path_p);
					if (path_p2 >= 0.1 && count2 > count) fprintf (stdout, ",AF=%f", path_p2);	
					fprintf (stdout, "\n");
				}
			}
		}
	}
}

