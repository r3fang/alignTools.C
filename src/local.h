/*--------------------------------------------------------------------*/
/* local.c                                                            */
/* Author: Rongxin Fang                                               */
/* Contact: r3fang@ucsd.edu                                           */
/* Date: 07-22-2015                                                   */
/* Pair wise local alignment without affine gap.                      */
/* initilize S(i,j):                                                  */
/* S(i, 0) = 0 and S(0, j) = 0, the first row and column to be 0      */
/* reccurrance relations:                                             */
/* S(i,j) = max{S(i-1, j-1)+s(x,y), S(i-1, j)+gap, S(i, j-1)+gap, 0}  */
/*               (S(i-1, j-1)+s(x,y))   # DIAGONAL                    */
/* S(i,j) = max  (   S(i-1, j)+gap  )   # RIGHT                       */
/*               (   S(i, j-1)+gap  )   # LEFT                        */
/*               (        0         )   # HOME                        */
/* Traceback:                                                         */
/* find the max value of S, can be anywhere in the matrix, start      */
/* tracing back from there and stop when we get to a cell = 0         */
/*--------------------------------------------------------------------*/
#ifndef _LOCAL_
#define _LOCAL_

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
trace_back_local(matrix_t *S, kstring_t *ks1, kstring_t *ks2, kstring_t *res_ks1, kstring_t *res_ks2, int i, int j){
	if(S == NULL || ks1 == NULL || ks2 == NULL || res_ks1 == NULL || res_ks2 == NULL) die("trace_back: parameter error");
	int m = 0; 
	while(S->M[i][j] != 0){
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
			case HOME: // go home
				i = 0;
				j = 0;
				break;
			default:
				break;	
		}
	}
	res_ks1->l = m;
	res_ks2->l = m;	
	/* reverse the string */
	res_ks1->s = strrev(res_ks1->s);
	res_ks2->s = strrev(res_ks2->s);
}
/*
 * main function for alignment	
 */
static inline double 
align_local(kstring_t *s1, kstring_t *s2, kstring_t *r1, kstring_t *r2){
	if(s1 == NULL || s2 == NULL || r1 == NULL || r2 == NULL) die("global: parameter error\n");
	size_t m   = s1->l + 1;
	size_t n   = s2->l + 1;
	matrix_t *S = create_matrix(m, n);
	size_t i, j, k, l;
	double max_score = -INFINITY;
	int i_max, j_max;
	int idx;
	for(i = 1; i <= s1->l; i++){
		for(j = 1; j <= s2->l; j++){
			double new_score = (strncmp(s1->s+(i-1), s2->s+(j-1), 1) == 0) ? MATCH : MISMATCH;
			idx = max5(&S->M[i][j], S->M[i][j-1] + GAP, S->M[i-1][j-1] + new_score, S->M[i-1][j] + GAP, 0.0, -INFINITY);
			if(idx == 0) S->pointerM[i][j] = LEFT;
			if(idx == 1) S->pointerM[i][j] = DIAGONAL;
			if(idx == 2) S->pointerM[i][j] = RIGHT;
			if(idx == 3) S->pointerM[i][j] = HOME;		
			if(max_score < S->M[i][j]){
				max_score = S->M[i][j];
				i_max = i; j_max = j;
			}
		}
	}
	// find max value of S->score, can be anywhere in matrix
	// stop when we get a cell with 0
	trace_back_local(S, s1, s2, r1, r2, i_max, j_max);
	destory_matrix(S);
	return max_score;
}

/* main function. */
static inline int
main_local(int argc, char *argv[]) {
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
	printf("%f\n", align_local(ks1, ks2, r1, r2));
	printf("%s\n%s\n", r1->s, r2->s);
	kstring_destory(ks1);
	kstring_destory(ks2);
	kstring_destory(r1);
	kstring_destory(r2);
	return 0;
}
#endif