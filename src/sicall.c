/*
 * sicall.c: Calculate the likelihoods of SNPs and INDELs. Make the calling.
 * Author: Mengyao Zhao
 * Create date: 2011-08-09
 * Contact: zhangmp@bc.edu
 * Last revise: 2012-06-05 
 */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "sicall.h"

void likelihood (double** transition, double** emission, char* ref, char* target_name, int32_t begin, int32_t filter)
{
	int32_t k, ref_len = strlen (ref);
//	for (k = 6; k < ref_len - 4; k ++) {	// for sliding window
	for (k = 1; k < ref_len + 1; ++k) {	// for small example test
		if (ref[k - 1] == 'A' || ref[k - 1] == 'a' || ref[k - 1] == 'C' || ref[k - 1] == 'c' || ref[k - 1] == 'G' || 
		ref[k - 1] == 'g' || ref[k - 1] == 'T' || ref[k - 1] == 't'/* || ref[k - 1] == 'N' || ref[k - 1] == 'n'*/) {

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
				default:
					fprintf(stderr, "Wrong reference sequence. \n");
					exit (1);
					break;
			}
			
			if (transition[k - 1][0] >= 0.1 && e <= 0.9) {
			//	double t = transition[k][0];
				float qual = -4.343 * log(e);
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
				if (base == 'N') goto indel;
				
				if (max > 0.9) {	// non reference allele
			//		qual = -4.343 * log (max);
					fprintf (stdout, "%s\t", target_name);
					fprintf (stdout, "%d\t.\t%c\t", begin + k, ref[k - 1]);
					fprintf (stdout, "%c\t%f\t", base, qual);
					if (filter == 0) fprintf (stdout, ".\t");
					else if (qual >= filter)	fprintf (stdout, "PASS\t");
					else fprintf (stdout, "q%d\t", filter);
					if (k == 1) fprintf (stdout, "AF=%f\n", max);
					else fprintf (stdout, "AF=%f\n", max * transition[k - 1][0]);
				} else {
					int32_t num2, n1, n2;
					double max2, temp1 = 0, temp2 = 0;
					char base2, b1 = 'N', b2 = 'N';
					if (num == 1) {
						if (emission[k][2] > emission[k][4]) {
							temp1 = emission[k][2];
							n1 = 2;
							b1 = 'C';
						} else {
							temp1 = emission[k][4];
							n1 = 4;
							b1 = 'G';
						}
					} else if (num == 2) {
						if (emission[k][1] > emission[k][4]) {
							temp1 = emission[k][1];
							n1 = 1;
							b1 = 'A';
						} else {
							temp1 = emission[k][4];
							n1 = 4;
							b1 = 'G';
						}
					} else if (num == 4 || num == 8) {
						if (emission[k][1] > emission[k][2]) {
							temp1 = emission[k][1];
							n1 = 1;
							b1 = 'A';
						} else {
							temp1 = emission[k][2];
							n1 = 2;
							b1 = 'C';
						}
					}
					if (num == 1 || num == 2 || num == 4) {
						if (emission[k][8] > emission[k][15]) {
							temp2 = emission[k][8];
							n2 = 8;
							b2 = 'T';
						} else {
							temp2 = emission[k][15];
							n2 = 15;
							b2 = 'N';
						}
					} else if (num == 8) {
						if (emission[k][4] > emission[k][15]) {
							temp2 = emission[k][4];
							n2 = 4;
							b2 = 'G';
						} else {
							temp2 = emission[k][15];
							n2 = 15;
							b2 = 'N';
						}
					}
					if (temp1 > temp2) {
						max2 = temp1;
						num2 = n1;
						base2 = b1;
					} else {
						max2 = temp2;
						num2 = n2;
						base2 = b2;
					}
					if (base != ref[k - 1]) {
						fprintf (stdout, "%s\t", target_name);
						fprintf (stdout, "%d\t.\t%c\t", begin + k, ref[k - 1]);
						if (max2 >= 0.1 && base2 != ref[k - 1] && base2 != 'N') fprintf (stdout, "%c,%c\t", base, base2);
						else fprintf (stdout, "%c\t", base);
						fprintf (stdout, "%f\t", qual);
						if (filter == 0) fprintf (stdout, ".\t");
						else if (qual >= filter)	fprintf (stdout, "PASS\t");
						else fprintf (stdout, "q%d\t", filter);
						if (max2 >= 0.1 && base2 != ref[k - 1] && base2 != 'N') fprintf (stdout, "AF=%f,AF=%f\n", max * transition[k - 1][0], max2 * transition[k - 1][0]);
						else fprintf (stdout, "AF=%f\n", max * transition[k - 1][0]);
					} else if (base2 != 'N') {
						fprintf (stdout, "%s\t", target_name);
						fprintf (stdout, "%d\t.\t%c\t", begin + k, ref[k - 1]);
						fprintf (stdout, "%c", base2);
						char base3, b;
						double max3, m;
						if ((max + max2) <= 0.9) {
							if (num + num2 == 3) {
								if (emission[k][4] > emission[k][8]) {
									m = emission[k][4];
									b = 'G';
								} else {
									m = emission[k][8];
									b = 'T';
								}
							} else if (num + num2 == 5) {
								if (emission[k][2] > emission[k][8]) {
									m = emission[k][2];
									b = 'C';
								} else {
									m = emission[k][8];
									b = 'T';
								}
							} else if (num + num2 == 9) {
								if (emission[k][2] > emission[k][4]) {
									m = emission[k][2];
									b = 'C';
								} else {
									m = emission[k][4];
									b = 'G';
								}
							} else if (num + num2 == 6) {
								if (emission[k][1] > emission[k][8]) {
									m = emission[k][1];
									b = 'A';
								} else {
									m = emission[k][8];
									b = 'T';
								}
							} else if (num + num2 == 10) {
								if (emission[k][1] > emission[k][4]) {
									m = emission[k][1];
									b = 'A';
								} else {
									m = emission[k][4];
									b = 'G';
								}
							} else if (num + num2 == 12) {
								if (emission[k][1] > emission[k][2]) {
									m = emission[k][1];
									b = 'A';
								} else {
									m = emission[k][2];
									b = 'C';
								}
							}
							if (m > emission[k][15]) {
								max3 = emission[k][15];
								base3 = 'N';
							}
							if (max3 >= 0.1 && base3 != 'N') fprintf (stdout, ",%c\t", base3);
							else fprintf (stdout, "\t");
						} else fprintf (stdout, "\t");	
						fprintf (stdout, "%f\t", qual);
						if (filter == 0) fprintf (stdout, ".\t");
						else if (qual >= filter)	fprintf (stdout, "PASS\t");
						else fprintf (stdout, "q%d\t", filter);
						if (max3 >= 0.1 && base3 != 'N' && k == 1) fprintf (stdout, "AF=%f,AF=%f\n", max2, max3);
						else if (max3 >= 0.1 && base3 != 'N' && k > 1) 
							fprintf (stdout, "AF=%f,AF=%f\n", max2 * transition[k - 1][0], max3 * transition[k - 1][0]);
						else if ((max3 < 0.1 || base3 == 'N') && k == 1) fprintf (stdout, "AF=%f\n", max2);
						else if ((max3 < 0.1 || base3 == 'N') && k > 1) fprintf (stdout, "AF=%f\n", max2 * transition[k - 1][0]);
					}
				}
			}

indel:
			if (transition[k][0] <= 0.9) {
				
				/* Detect insertion. */
				int32_t insert_num = 0;
				double insert_p = 0, p;
				if (emission[k][0] >= 0.1) {
					insert_p += emission[k][0];
					insert_num ++;	//FIXME: the insert num calculation is wrong.
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

				p = transition[k][1] * (transition[k][4] + transition[k][5]);

				if (p >= 0.1) {
					double t = transition[k][0];
					float qual = -4.343 * log(t/(t + p));
					fprintf (stdout, "%s\t%d\t.\t%c\t<I%d>\t%f\t", target_name, begin + k, ref[k - 1], insert_num, qual);
					if (filter == 0) fprintf (stdout, ".\t");
					else if (qual >= filter)	fprintf (stdout, "PASS\t");
					else fprintf (stdout, "q%d\t", filter);
					fprintf (stdout, "AF=%f\n", p);	//FIXME: AF calculation is wrong.
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
						int32_t count_p = path_p2 > path_p ? count2 : count;
						double p = path_p2 > path_p ? path_p2 : path_p;
						for (i = 1; i <= count_p; ++i) t *= transition[k + i][0];
				//		t /= (t + p);
						float qual = -4.343 * log(t/(t + p));
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

