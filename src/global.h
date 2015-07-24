/*--------------------------------------------------------------------*/
/* global.h 		                                                  */
/* Author: Rongxin Fang                                               */
/* E-mail: r3fang@ucsd.edu                                            */
/* Date: 07-22-2015                                                   */
/* Pair wise global alignment without affine gap.                     */
/* initilize S(i,j):                                                  */
/* S(i, 0) = g*i; S(0, j) = g*j                                       */
/* reccurrance relations:                                             */
/* S(i,j) = max{S(i-1, j-1)+s(x,y), S(i-1, j)+gap, S(i, j-1)+gap}     */
/*               (S(i-1, j-1) +s(x,y))   # DIAGONAL                   */
/* S(i,j) = max  (   S(i-1, j) + gap  )   # RIGHT                     */
/*               (   S(i, j-1) + gap  )   # LEFT                      */
/* Traceback:                                                         */
/* S(m, n) holds the optimal alignment score; trace pointers back     */
/* from S(m, n) to S(0, 0) to recover alignment.                      */
/*--------------------------------------------------------------------*/
#ifndef _GLOBAL_
#define _GLOBAL_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "utils.h"     // die and mycalloc
#include "kseq.h"      // fasta parser
#include "kstring.h"   // kstring_t type
#include "alignment.h" // basic functions

static inline void 
trace_back_gl(matrix_t *S, kstring_t *ks1, kstring_t *ks2, kstring_t *res_ks1, kstring_t *res_ks2, int i, int j){
	if(S == NULL || ks1 == NULL || ks2 == NULL || res_ks1 == NULL || res_ks2 == NULL) die("trace_back: parameter error");
	int m = 0; 
	while(i>0 && j >0){
		switch(S->pointerM[i][j]){
			case LEFT:
				res_ks2->s[m] = ks2->s[--j];
				res_ks1->s[m++] = '-';
				break;
			case DIAGONAL:
				res_ks1->s[m] = ks1->s[--i];
				res_ks2->s[m++] = ks2->s[--j];
				break;
			case RIGHT:
				res_ks1->s[m] = ks1->s[--i];
				res_ks2->s[m++] = '-';
				break;
			default:
				break;	
		}
	}
	if(j>0){while(j>0){
			res_ks1->s[m] = '-';	
			res_ks2->s[m++] = ks2->s[--j];				
		}
	}
	if(i>0){while(i>0){
			res_ks2->s[m] = '-';	
			res_ks1->s[m++] = ks1->s[--i];				
		}
	}	
	res_ks1->l = m;
	res_ks2->l = m;	
	res_ks1->s = strrev(res_ks1->s);
	res_ks2->s = strrev(res_ks2->s);
}

static inline double 
align_gl(kstring_t *s1, kstring_t *s2, kstring_t *r1, kstring_t *r2){
	if(s1 == NULL || s2 == NULL || r1 == NULL || r2 == NULL) die("global: parameter error\n");
	size_t m   = s1->l + 1;
	size_t n   = s2->l + 1;
	matrix_t *S = create_matrix(m, n);
	// initlize the first row and column of S
	size_t i, j, k, l;
	for(i=0; i < S->m; i++) S->M[i][0] = GAP*i;
	for(j=0; j < S->n; j++) S->M[0][j] = GAP*j;
	
	int idx;
	for(i = 1; i <= s1->l; i++){
		for(j = 1; j <= s2->l; j++){
			double new_score = (strncmp(s1->s+(i-1), s2->s+(j-1), 1) == 0) ? MATCH : MISMATCH;
			idx= max5(&S->M[i][j], S->M[i][j-1] + GAP, S->M[i-1][j-1] + new_score, S->M[i-1][j] + GAP, -INFINITY, -INFINITY);
			if(idx==0) S->pointerM[i][j] = LEFT;
			if(idx==1) S->pointerM[i][j] = DIAGONAL;
			if(idx==2) S->pointerM[i][j] = RIGHT;
		}
	}
	trace_back_gl(S, s1, s2, r1, r2, s1->l, s2->l);
	double res = S->M[s1->l][s2->l];
	destory_matrix(S);
	return res;
}

/* main function. */
static inline int 
main_global(int argc, char *argv[]) {
	kstring_t *ks1, *ks2; 
	ks1 = mycalloc(1, kstring_t);
	ks2 = mycalloc(1, kstring_t);
	if (argc == 1) {
		fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
		return 1;
	}
	kstring_read(argv[1], ks1, ks2);
	if(ks1->s == NULL || ks2->s == NULL) die("fail to read sequence\n");
	kstring_t *r1 = mycalloc(1, kstring_t);
	kstring_t *r2 = mycalloc(1, kstring_t);
	r1->s = mycalloc(ks1->l + ks2->l, char);
	r2->s = mycalloc(ks1->l + ks2->l, char);
	printf("score=%f\n", align_gl(ks1, ks2, r1, r2));
	printf("%s\n%s\n", r1->s, r2->s);
	kstring_destory(ks1);
	kstring_destory(ks2);
	kstring_destory(r1);
	kstring_destory(r2);
	return 0;
}
#endif