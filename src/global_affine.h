/*--------------------------------------------------------------------*/
/* global_affine_gap.h 	                                              */
/* Author: Rongxin Fang                                               */
/* E-mail: r3fang@ucsd.edu                                            */
/* Date: 07-22-2015                                                   */
/* Pair wise global alignment with affine gap.                        */
/* initilize L(i,j), M(i,j), U(i,j):                                  */
/* M(0,0) = 0                                                         */
/* M(i,0) = -INF               (i>0 )                                 */
/* M(0,j) = -INF               (j>0 )                                 */
/* L(i,0) = open + extension*i (i>=0)                                 */
/* L(0,j) = -INF               (j>0 )                                 */
/* U(0,j) = open + extension*j (j>=0)                                 */
/* reccurrance relations:                                             */
/* M(i,j) = max{M(i-1, j-1), U(i-1, j-1), L(i-1, j-1)}+s(x,y)         */
/* U(i,j) = max{M(i-1, j)+gap, U(i-1, j)+extension}                   */
/* L(i,j) = max{M(i, j-1)+gap, L(i, j-1)+extension}                   */
/* Traceback:                                                         */
/* start at largest of M(m,n), L(m,n), U(m,n)                         */
/* Stop at any of M(0,0), I(0,0), U(0,0)                              */
/*--------------------------------------------------------------------*/
#ifndef _GLOBAL_AFFINE_
#define _GLOBAL_AFFINE_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "utils.h"
#include "kseq.h"
#include "kstring.h"
#include "alignment.h"

static inline void 
trace_back_gla(matrix_t *S, kstring_t *s1, kstring_t *s2, kstring_t *res_ks1, kstring_t *res_ks2, int state){
	if(S == NULL || s1 == NULL || s2 == NULL || res_ks1 == NULL || res_ks2 == NULL) die("trace_back: paramter error");
	int i = s1->l; int j = s2->l;
	int cur = 0; 
	while(i > 0 && j > 0){
		switch(state){
			case LOW:
				state = S->pointerL[i][j]; // change to next state
				res_ks1->s[cur] = s1->s[--i];
				res_ks2->s[cur++] = '-';
				break;
			case MID:
				state = S->pointerM[i][j]; // change to next state
                res_ks1->s[cur] = s1->s[--i];
                res_ks2->s[cur++] = s2->s[--j];
				break;
			case UPP:
				state = S->pointerU[i][j];
				res_ks1->s[cur] = '-';
            	res_ks2->s[cur++] = s2->s[--j];
				break;
			default:
				break;
			}
	}	
	if(j>0){while(j>0){
			res_ks1->s[cur] = '-';	
			res_ks2->s[cur++] = s2->s[--j];				
		}
	}
	if(i>0){while(i>0){
			res_ks2->s[cur] = '-';	
			res_ks1->s[cur++] = s1->s[--i];				
		}
	}	
	res_ks1->l = cur;
	res_ks2->l = cur;	
	res_ks1->s = strrev(res_ks1->s);
	res_ks2->s = strrev(res_ks2->s);
}

/*
 * Global alignment with affine gap penality
 */
static inline double 
align_gla(kstring_t *s1, kstring_t *s2, kstring_t *r1, kstring_t *r2){
	if(s1 == NULL || s2 == NULL || r1 == NULL || r2 == NULL) die("align: parameter error\n");
	size_t m   = s1->l + 1; size_t n   = s2->l + 1;
	matrix_t *S = create_matrix(m, n);
	// initlize DP matrix
	S->M[0][0] = 0.0;
	S->L[0][0] = S->U[0][0] = GAP;
	// initlize 0 column
	int i, j;
	for(i=1; i<S->m; i++){
		S->L[i][0] = GAP + EXTENSION*(i);
		S->M[i][0] = -INFINITY;
		S->U[i][0] = -INFINITY;
	}
	for(j=1; j<S->n; j++){
		S->L[0][j] = -INFINITY;
		S->M[0][j] = -INFINITY;
		S->U[0][j] = GAP + EXTENSION*(j);
	}
	//-------------------------------
	double new_score;
	int idx;
	// recurrance relation
	for(i=1; i<=s1->l; i++){
		for(j=1; j<=s2->l; j++){
			// MID
			new_score = (strncmp(s1->s+(i-1), s2->s+(j-1), 1) == 0) ? MATCH : MISMATCH;
			//new_score = match(s1->s[i-1], s2->s[j-1], BLOSUM62);
			idx = max5(&S->M[i][j], S->L[i-1][j-1]+new_score, S->M[i-1][j-1]+new_score, S->U[i-1][j-1]+new_score, -INFINITY, -INFINITY);
			if(idx==0) S->pointerM[i][j] = LOW;
			if(idx==1) S->pointerM[i][j] = MID;
			if(idx==2) S->pointerM[i][j] = UPP;
			// LOW
			idx = max5(&S->L[i][j], S->L[i-1][j]+EXTENSION, S->M[i-1][j]+GAP, -INFINITY, -INFINITY, -INFINITY);
			if(idx==0) S->pointerL[i][j] = LOW;
			if(idx==1) S->pointerL[i][j] = MID;			
			// UPP
			idx = max5(&S->U[i][j], -INFINITY, S->M[i][j-1]+GAP, S->U[i][j-1]+EXTENSION, -INFINITY, -INFINITY);
			if(idx==1) S->pointerU[i][j] = MID;
			if(idx==2) S->pointerU[i][j] = UPP;
		}
	}
	double max_score; int max_state;
	idx = max5(&max_score, S->L[s1->l][s2->l], S->M[s1->l][s2->l], S->U[s1->l][s2->l], -INFINITY, -INFINITY);
	if(idx==0) max_state = LOW;
	if(idx==1) max_state = MID;
	if(idx==2) max_state = UPP;
	trace_back_gla(S, s1, s2, r1, r2, max_state);	
	destory_matrix(S);
	return max_score;
}

/* main function. */
static inline int 
main_global_affine(int argc, char *argv[]) {
	opt_t *opt = NULL;	
	kstring_t *ks1, *ks2; 
	ks1 = mycalloc(1, kstring_t);
	ks2 = mycalloc(1, kstring_t);
	if (argc == 1) {
		fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
		return 1;
	}
	kstring_read(argv[1], ks1, ks2, opt);
	if(ks1->s == NULL || ks2->s == NULL) die("fail to read sequence\n");
	kstring_t *r1 = mycalloc(1, kstring_t);
	kstring_t *r2 = mycalloc(1, kstring_t);
	r1->s = mycalloc(ks1->l + ks2->l, char);
	r2->s = mycalloc(ks1->l + ks2->l, char);
	printf("score=%f\n", align_gla(ks1, ks2, r1, r2));
	printf("%s\n%s\n", r1->s, r2->s);
	kstring_destory(ks1);
	kstring_destory(ks2);
	kstring_destory(r1);
	kstring_destory(r2);
	return 0;
}
#endif