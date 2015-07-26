/*--------------------------------------------------------------------*/
/* local_affine.h	                                                  */
/* Author: Rongxin Fang                                               */
/* E-mail: r3fang@ucsd.edu                                            */
/* Date: 07-22-2015                                                   */
/* Pair wise local alignment with affine gap.                         */
/* initilize L(i,j), M(i,j), U(i,j):                                  */
/* M(0,0) = 0                                                         */
/* M(i,0) = 0                  (i>0 )                                 */
/* M(0,j) = 0                  (j>0 )                                 */
/* L(0,j) = -INF               (j>=0)                                 */
/* U(i,0) = -INF               (j>=0)                                 */
/* reccurrance relations:                                             */
/* M(i,j) = max{M(i-1, j-1)+s(x,y), U(i-1, j-1)+s(x,y),               */
/*              L(i-1, j-1)+s(x,y), 0}                                */
/* U(i,j) = max{M(i-1, j)+gap, U(i-1, j)+extension}                   */
/* L(i,j) = max{M(i, j-1)+gap, L(i, j-1)+extension}                   */
/* Traceback:                                                         */
/* start at max(M(i,j))                                               */
/* Stop at any of M(0,0)                                              */
/*--------------------------------------------------------------------*/
#ifndef _LOCAL_AFFINE_
#define _LOCAL_AFFINE_

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
trace_back_local_affine(matrix_t *S, kstring_t *s1, kstring_t *s2, kstring_t *res_ks1, kstring_t *res_ks2, int i, int j){
	if(S == NULL || s1 == NULL || s2 == NULL || res_ks1 == NULL || res_ks2 == NULL) die("trace_back: paramter error");
	int cur = 0; 
	int state = MID;
	while(i>0 && j>0){
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
			case HOME:
				i = 0;
				j = 0;
				break;
			default:
				break;
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
align_local_affine(kstring_t *s1, kstring_t *s2, kstring_t *r1, kstring_t *r2){
	if(s1 == NULL || s2 == NULL || r1 == NULL || r2 == NULL) die("align: parameter error\n");
	size_t m   = s1->l + 1; size_t n   = s2->l + 1;
	matrix_t *S = create_matrix(m, n);
	int i, j;
	int i_max, j_max;
	double max_score = -INFINITY;
	double new_score;
	int idx;
	// recurrance relation
	for(i=1; i<=s1->l; i++){
		for(j=1; j<=s2->l; j++){
			// MID
			new_score = ((s1->s[i-1] - s2->s[j-1]) == 0) ? MATCH : MISMATCH;
			idx = max5(&S->M[i][j], S->L[i-1][j-1]+new_score, S->M[i-1][j-1]+new_score, S->U[i-1][j-1]+new_score, 0.0, -INFINITY);
			if(idx==0) S->pointerM[i][j] = LOW;
			if(idx==1) S->pointerM[i][j] = MID;
			if(idx==2) S->pointerM[i][j] = UPP;
			if(idx==3) S->pointerM[i][j] = HOME;			
			if(S->M[i][j] > max_score){
				max_score = S->M[i][j];
				i_max = i; j_max = j;
			}
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
	trace_back_local_affine(S, s1, s2, r1, r2, i_max, j_max);	
	destory_matrix(S);
	return max_score;
}

/* main function. */
static inline int 
main_local_affine(int argc, char *argv[]) {
	opt_t *opt;
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
	printf("score=%f\n", align_local_affine(ks1, ks2, r1, r2));
	printf("%s\n%s\n", r1->s, r2->s);
	kstring_destory(ks1);
	kstring_destory(ks2);
	kstring_destory(r1);
	kstring_destory(r2);
	return 0;
}
#endif