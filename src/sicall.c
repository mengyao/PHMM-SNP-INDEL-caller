/*
 * sicall.c: Calculate the likelihoods of SNPs and INDELs. Make the calling.
 * Author: Mengyao Zhao
 * Create date: 2011-08-09
 * Contact: zhangmp@bc.edu
 * Last revise: 2011-10-03 
 */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "sicall.h"

void likelihood (double** transition, double** emission, char* ref, char* target_name, int32_t begin, int32_t filter)
{
	int32_t k, ref_len = strlen (ref);
	for (k = 3; k < ref_len - 1; k ++) {
		if (ref[k - 1] == 'A' || ref[k - 1] == 'a' || ref[k - 1] == 'C' || ref[k - 1] == 'c' || ref[k - 1] == 'G' || 
		ref[k - 1] == 'g' || ref[k - 1] == 'T' || ref[k - 1] == 't' || ref[k - 1] == 'N' || ref[k - 1] == 'n') {

			/* Detect SNP. */
			double e;
			switch (ref[k - 1]) {
				case 'A':
				case 'a':
					e = emission[k][1];
					break;
				case 'C':
				case 'c':
					e = emission[k][2];
					break;
				case 'G':
				case 'g':
					e = emission[k][4];
					break;
				case 'T':
				case 't':
					e = emission[k][8];
					break;
				case 'N':
				case 'n':
					e = emission[k][15];
					break;
				default:
					fprintf(stderr, "Wrong reference sequence. \n");
					exit (1);
					break;
			}
			
			if (transition[k - 1][0] >= 0.1 && e <= 0.9) {
				/*double t = k == ref_len ? (1 - transition[k][1]) : transition[k][0];*/
				double t = transition[k][0];
				float qual = -4.343 * log(e * t);
				double max;
				char base;
				int32_t num;
				if (emission[k][1] > emission[k][2]) {
					max = emission[k][1];
					base = 'A';
					num = 1;
				} else {
					max = emission[k][2];
					base = 'C';
					num = 2;
				}
				if (max < emission[k][4]) {
					max = emission[k][4];
					base = 'G';
					num = 4;
				}
				if (max < emission[k][8]) {
					max = emission[k][8];
					base = 'T';
					num = 8;
				}
				if (max < emission[k][15]) {
					max = emission[k][15];
					base = 'N';
					num = 15;
				}
				fprintf (stdout, "%s\t", target_name);
				fprintf (stdout, "%d\t.\t%c\t", begin + k, ref[k - 1]);
				if (max > 0.9) {
					fprintf (stdout, "%c\t%f\t", base, qual);
					if (filter == 0) fprintf (stdout, ".\t");
					else if (qual >= filter)	fprintf (stdout, "PASS\t");
					else fprintf (stdout, "q%d\t", filter);
					if (k == 1) fprintf (stdout, "AF=%f\n", max);
					else fprintf (stdout, "AF=%f\n", max * transition[k - 1][0]);
				} else if (max <= 0.9 && base != ref[k - 1]) {
					double max2 = 0;
					char base2 = 'N';
					if (emission[k][1] > emission[k][2] && num != 1) {
						max2 = emission[k][1];
						base2 = 'A';
					} else if (emission[k][1] <= emission[k][2] && num != 2) {
						max2 = emission[k][2];
						base2 = 'C';
					}
					if (max2 < emission[k][4] && num != 4) {
						max2 = emission[k][4];
						base2 = 'G';
					}
					if (max2 < emission[k][8] && num != 8) {
						max2 = emission[k][8];
						base2 = 'T';
					}
					if (max2 >= 0.1 && base2 != ref[k - 1]) fprintf (stdout, "%c,%c\t", base, base2);
					else fprintf (stdout, "%c\t", base);
					fprintf (stdout, "%f\t", qual);
					if (filter == 0) fprintf (stdout, ".\t");
					else if (qual >= filter)	fprintf (stdout, "PASS\t");
					else fprintf (stdout, "q%d\t", filter);
					if (max2 >= 0.1 && base2 != ref[k - 1]) fprintf (stdout, "AF=%f,AF=%f\n", max * transition[k - 1][0], max2 * transition[k - 1][0]);
					else fprintf (stdout, "AF=%f\n", max * transition[k - 1][0]);
				} else {
					double max2 = 0, max3 = 0;
					char base2 = 'N';
					int32_t num2 = 15;		
					if (emission[k][1] > emission[k][2] && num != 1) {
						max2 = emission[k][1];
						base2 = 'A';
						num2 = 1;
					} else if (emission[k][1] <= emission[k][2] && num != 2) {
						max2 = emission[k][2];
						base2 = 'C';
						num2 = 2;
					}
					if (max2 < emission[k][4] && num != 4) {
						max2 = emission[k][4];
						base2 = 'G';
						num2 = 4;
					}
					if (max2 < emission[k][8] && num != 8) {
						max2 = emission[k][8];
						base2 = 'T';
						num2 = 8;
					}
					fprintf (stdout, "%c", base2);
					if ((max + max2) <= 0.9) {
						char base3 = 'N';
						int32_t num3;		
						if (emission[k][1] > emission[k][2] && num != 1 && num2 != 1) {
							max3 = emission[k][1];
							base3 = 'A';
							num3 = 1;
						} else if (emission[k][1] <= emission[k][2] && num != 2 && num2 != 2) {
							max3 = emission[k][2];
							base3 = 'C';
							num3 = 2;
						}
						if (max3 < emission[k][4] && num != 4 && num2 != 4) {
							max3 = emission[k][4];
							base3 = 'G';
							num3 = 4;
						}
						if (max3 < emission[k][8] && num != 8 && num2 != 8) {
							max3 = emission[k][8];
							base3 = 'T';
							num3 = 8;
						}
						if (max3 >= 0.1) fprintf (stdout, ",%c\t", base3);
						else fprintf (stdout, "\t");
					} else fprintf (stdout, "\t");	
					fprintf (stdout, "%f\t", qual);
					if (filter == 0) fprintf (stdout, ".\t");
					else if (qual >= filter)	fprintf (stdout, "PASS\t");
					else fprintf (stdout, "q%d\t", filter);
					if (max3 >= 0.1 && k == 1) fprintf (stdout, "AF=%f,AF=%f\n", max2, max3);
					else if (max3 >= 0.1 && k > 1) 
						fprintf (stdout, "AF=%f,AF=%f\n", max2 * transition[k - 1][0], max3 * transition[k - 1][0]);
					else if (max3 < 0.1 && k == 1) fprintf (stdout, "AF=%f\n", max2);
					else fprintf (stdout, "AF=%f\n", max2 * transition[k - 1][0]);
				}
			}

			if (transition[k][0] <= 0.9) {
				
				/* Detect insertion. */
				int32_t insert_num = 0;
				double insert_p = 0, transition_p, p;
				if (emission[k][0] >= 0.1) {
					insert_p += emission[k][0];
					insert_num ++;
				}
				if (emission[k][3] >= 0.1) {
					insert_p += emission[k][3];
					insert_num ++;
				}
				if (emission[k][5] >= 0.1) {
					insert_p += emission[k][5];
					insert_num ++;
				}
				if (emission[k][9] >= 0.1) {
					insert_p += emission[k][9];
					insert_num ++;
				}
				if (emission[k][14] >= 0.1) {
					insert_p += emission[k][14];
					insert_num ++;
				}

				transition_p = transition[ref_len][5];
				p = transition[k][1] * insert_p * transition_p;

				if (p >= 0.1) {
					double t = transition[k][0];
					float qual = -4.343 * log(e * t);
					fprintf (stdout, "%s\t%d\t.\t%c\t<I%d>\t%f\t", target_name, begin + k, ref[k - 1], insert_num, qual);
					if (filter == 0) fprintf (stdout, ".\t");
					else if (qual >= filter)	fprintf (stdout, "PASS\t");
					else fprintf (stdout, "q%d\t", filter);
					fprintf (stdout, "AF=%f\n", p);
				}

				/* Detect deletion. */
				if (transition[k][2] >= 0.1) {
					double path_p = transition[k][2];
					double path_p2 = transition[k][2];
					int32_t count = 1, count2 = 1, flag = 0;
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
						double t = transition[k][0];
						float qual = -4.343 * log(e * t);
						fprintf (stdout, "%s\t%d\t.\t%c", target_name, begin + k, ref[k - 1]);
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
}

