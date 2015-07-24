/*--------------------------------------------------------------------*/
/* fit_affine_jump.h 	                                              */
/* Author: Rongxin Fang                                               */
/* E-mail: r3fang@ucsd.edu                                            */
/* Date: 07-23-2015                                                   */
/* Pair wise fit alignment with affine gap and jump state.            */
/* This could be used to align RNA-seq reads with intron splicing as  */
/* jump state.                                                        */
/*                                                                    */
/* initilize L(i,j), M(i,j), U(i,j):                                  */
/*--------------------------------------------------------------------*/
/* M(0,0) = -INF                                                      */
/* M(i,0) = -INF               (i>0 )                                 */
/* M(0,j) = 0                  (j>0 )                                 */
/* L(0,j) = -INF               (j>=0)                                 */
/* L(i,0) = -INF               (i>=0)                                 */
/* U(i,0) = -INF               (i>0 )                                 */
/* U(0,j) = 0                  (j>=0)                                 */
/* J(0,j) = -INF               (j>=0)                                 */
/* J(i,0) = -INF               (i>=0)                                 */
/*                                                                    */
/* reccurrance relations:                                             */
/*--------------------------------------------------------------------*/
/* M(i,j) = max{M(i-1, j-1)+s(x,y), U(i-1, j-1)+s(x,y),               */
/*              L(i-1, j-1)+s(x,y), J(i-1, j-1)+s(x,y)}               */
/* U(i,j) = max{M(i-1, j)+GAP, U(i-1, j)+EXTENSION}                   */
/* L(i,j) = max{M(i, j-1)+GAP, L(i, j-1)+EXTENSION}                   */
/* We allow pointer move from M to J only at given positions on s2    */
/* J(i,j) = max{M(i-1, j)+JUMP, U(i-1, j)} if s2[j] = junction   OR   */
/* J(i,j) = max{M(i-1, j)-INFINITY, U(i-1, j)} if s2[j] = junction    */
/* Traceback:                                                         */
/*--------------------------------------------------------------------*/
/* start at max(M(m,j_max), L(n, j_max)), Stop at any of i=0 on M/L;  */
/* The rational behind this is no gap allowed to flank s1             */
/*--------------------------------------------------------------------*/
/* NOTE:                                                              */
/* 1. WE ONLY ALLOW STATE CHANGE FROM M(MATCH) TO J(JUMP) AT SPECIFIC */
/* POSITIONS ON S2.                                                   */ 
/* 2. WE ALLOW STATE CHANGE FROM JUMP TO MATCH AT ANYWHERE            */
/* 3. THE READ SHOULD NOT START WITH JUMP AND ALSO SHOULD NOT END     */
/* WITH JUMP STATE.                                                   */ 
/*--------------------------------------------------------------------*/
#ifndef _FIT_AFFINE_JUMP_
#define _FIT_AFFINE_JUMP_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include "utils.h"
#include "kseq.h"
#include "kstring.h"
#include "alignment.h"

//junctions
typedef struct {
	size_t size;
	int *pos;
} junction_t;

static inline void 
trace_back_fit_affine_jump(matrix_t *S, kstring_t *s1, kstring_t *s2, kstring_t *res_ks1, kstring_t *res_ks2, int state, int i, int j){
	if(S == NULL || s1 == NULL || s2 == NULL || res_ks1 == NULL || res_ks2 == NULL) die("trace_back: paramter error");
	int cur = 0; 
	while(i>0){
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
			case JUMP:
				state = S->pointerJ[i][j];
				res_ks1->s[cur] = '-';
	           	res_ks2->s[cur++] = s2->s[--j];
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
 * fit alignment with affine gap penality
 */
static inline double 
align_fit_affine_jump(kstring_t *s1, kstring_t *s2, kstring_t *r1, kstring_t *r2, junction_t *junctions){
	if(s1 == NULL || s2 == NULL || r1 == NULL || r2 == NULL) die("align: parameter error\n");
	if(s1->l > s2->l) die("first sequence must be shorter than the second to do fitting alignment"); 
	size_t m   = s1->l + 1; size_t n   = s2->l + 1;
	matrix_t *S = create_matrix(m, n);
	// initlize leftmost column
	int i, j;
	for(i=0; i<S->m; i++){
		S->M[i][0] = -INFINITY;
		S->U[i][0] = -INFINITY;
		S->L[i][0] = -INFINITY;
		S->J[i][0] = -INFINITY;
	}
	// initlize first row
	for(j=0; j<S->n; j++){
		S->M[0][j] = 0.0;
		S->U[0][j] = 0.0;
		S->J[0][j] = -INFINITY;
		S->L[0][j] = -INFINITY;
	}
	double new_score;
	int idx;
	
	// recurrance relation
	for(i=1; i<=s1->l; i++){
		for(j=1; j<=s2->l; j++){
			// MID any state can goto MID
			new_score = (strncmp(s1->s+(i-1), s2->s+(j-1), 1) == 0) ? MATCH : MISMATCH;
			idx = max5(&S->M[i][j], S->L[i-1][j-1]+new_score, S->M[i-1][j-1]+new_score, S->U[i-1][j-1]+new_score, S->J[i-1][j-1]+new_score, -INFINITY);
			if(idx == 0) S->pointerM[i][j]=LOW;
			if(idx == 1) S->pointerM[i][j]=MID;
			if(idx == 2) S->pointerM[i][j]=UPP;
			if(idx == 3) S->pointerM[i][j]=JUMP;			 
			
			// LOW
			idx = max5(&S->L[i][j], S->L[i-1][j]+EXTENSION, S->M[i-1][j]+GAP, -INFINITY, -INFINITY, -INFINITY);
			if(idx == 0) S->pointerL[i][j]=LOW;
			if(idx == 1) S->pointerL[i][j]=MID;
			
			// UPP
			idx = max5(&S->U[i][j], -INFINITY, S->M[i][j-1]+GAP, S->U[i][j-1]+EXTENSION, -INFINITY, -INFINITY);
			if(idx == 1) S->pointerU[i][j]=MID;
			if(idx == 2) S->pointerU[i][j]=UPP;
			
			// JUMP only allowed going to JUMP state at junction sites
			if(isvalueinarray(j-1, junctions->pos, junctions->size)){
				idx = max5(&S->J[i][j], -INFINITY, S->M[i][j-1]+JUMP_PENALITY, -INFINITY, S->J[i][j-1], -INFINITY);
				if(idx == 1) S->pointerJ[i][j] = MID;			
				if(idx == 3) S->pointerJ[i][j] = JUMP;			
			}else{
				idx = max5(&S->J[i][j], -INFINITY, -INFINITY, -INFINITY, S->J[i][j-1], -INFINITY);				
				if(idx == 3) S->pointerJ[i][j] = JUMP;
			}
		}
	}
	
	// find trace-back start point
	int i_max, j_max;
	double max_score = -INFINITY;
	int max_state;
	i_max = s1->l;
	for(j=0; j<s2->l; j++){
		if(max_score < S->M[i_max][j]){
			max_score = S->M[i_max][j];
			j_max = j;
			max_state = MID;
		}
	}
	for(j=0; j<s2->l; j++){
		if(max_score < S->L[i_max][j]){
			max_score = S->L[i_max][j];
			j_max = j;
			max_state = LOW;
		}
	}
	trace_back_fit_affine_jump(S, s1, s2, r1, r2, max_state, i_max, j_max);	
	destory_matrix(S);
	return max_score;
}

/* main function. */
static inline int 
main_fit_affine_jump(int argc, char *argv[]) {
	junction_t *junctions = mycalloc(1, junction_t);
	junctions->size = 5;
	junctions->pos = mycalloc(junctions->size, int);
	junctions->pos[0] = 90;  
	junctions->pos[1] = 91;  
	junctions->pos[2] = 92;  
	junctions->pos[4] = 89;  
	junctions->pos[5] = 558;
	junctions->pos[6] = 559;
	junctions->pos[7] = 560;
	
	kstring_t *ks1, *ks2; 
	ks1 = mycalloc(1, kstring_t);
	ks2 = mycalloc(1, kstring_t);
	if (argc == 1) {
		fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
		return 1;
	}
	kstring_read(argv[1], ks1, ks2);
	if(ks1->s == NULL || ks2->s == NULL) die("fail to read sequence\n");
	if(ks1->l > ks2->l) die("first sequence must be shorter than the second\n");
	kstring_t *r1 = mycalloc(1, kstring_t);
	kstring_t *r2 = mycalloc(1, kstring_t);
	r1->s = mycalloc(ks1->l + ks2->l, char);
	r2->s = mycalloc(ks1->l + ks2->l, char);
	printf("score=%f\n", align_fit_affine_jump(ks1, ks2, r1, r2, junctions));
	printf("%s\n%s\n", r1->s, r2->s);
	kstring_destory(ks1);
	kstring_destory(ks2);
	kstring_destory(r1);
	kstring_destory(r2);
	return 0;
}
#endif