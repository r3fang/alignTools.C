/*--------------------------------------------------------------------*/
/* global_affine_gap.c 	                                              */
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "utils.h"
#include "kseq.h"
#include "kstring.h"
#include "alignment.h"

//--------------
#define GL_ERR_NONE              0

// penality
//--------------
#define GAP                     -11.0
#define EXTENSION               -1.0
#define MATCH                    2.0
#define MISMATCH                -1.0

// state
//--------------
#define LOW                     100
#define MID                     200
#define UPP                     300

// DP matrix
typedef struct {
  unsigned int m;
  unsigned int n;
  double **L;
  double **M;
  double **U;
  int **pointerL;
  int **pointerM;
  int **pointerU;
} matrix_t;

/* store scoring matrix */
scoring_matrix_t *BLOSUM62;

/*
 * initilize matrix.	
 */
matrix_t *create_matrix(size_t m, size_t n){
	size_t i, j; 
	matrix_t *S = mycalloc(1, matrix_t);
	S->m = m;
	S->n = n;
	S->L = mycalloc(m, double*);
	S->M = mycalloc(m, double*);
	S->U = mycalloc(m, double*);
	S->pointerL = mycalloc(m, int*);
	S->pointerM = mycalloc(m, int*);
	S->pointerU = mycalloc(m, int*);
	for (i = 0; i < m; i++) {
		S->L[i] = mycalloc(n, double);
		S->M[i] = mycalloc(n, double);
		S->U[i] = mycalloc(n, double);
		S->pointerL[i] = mycalloc(n, int);
		S->pointerM[i] = mycalloc(n, int);
		S->pointerU[i] = mycalloc(n, int);
	}
	for(i=0; i<S->m; i++){
		for(j=0; j<S->n; j++){
			S->L[i][j] = -INFINITY;
			S->M[i][j] = -INFINITY;
			S->U[i][j] = -INFINITY;
		}
	}
	
	// initlize DP matrix
	S->M[0][0] = 0.0;
	S->L[0][0] = S->U[0][0] = GAP;
	// initlize 0 column
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
	return S;
}

/*
 * destory matrix
 */
int destory_matrix(matrix_t *S){
	if(S == NULL) die("destory_matrix: parameter error\n");
	int i, j;
	for(i = 0; i < S->m; i++){
		free(S->L[i]);
		free(S->M[i]);
		free(S->U[i]);
		free(S->pointerL[i]);
		free(S->pointerM[i]);
		free(S->pointerU[i]);
	}
	free(S);
	return GL_ERR_NONE;
}

int max3(double *res, double a1, double a2, double a3){
	*res = -INFINITY;
	int state;
	if(a1 > *res){*res = a1; state = LOW;}
	if(a2 > *res){*res = a2; state = MID;}
	if(a3 > *res){*res = a3; state = UPP;}
	return state;
}

void trace_back(matrix_t *S, kstring_t *s1, kstring_t *s2, kstring_t *res_ks1, kstring_t *res_ks2, int state){
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
double align(kstring_t *s1, kstring_t *s2, kstring_t *r1, kstring_t *r2){
	if(s1 == NULL || s2 == NULL || r1 == NULL || r2 == NULL) die("align: parameter error\n");
	size_t m   = s1->l + 1; size_t n   = s2->l + 1;
	matrix_t *S = create_matrix(m, n);
	int i, j;
	double new_score;
	// recurrance relation
	for(i=1; i<=s1->l; i++){
		for(j=1; j<=s2->l; j++){
			// MID
			//new_score = (strncmp(s1->s+(i-1), s2->s+(j-1), 1) == 0) ? MATCH : MISMATCH;
			new_score = match(s1->s[i-1], s2->s[j-1], BLOSUM62);
			S->pointerM[i][j] = max3(&S->M[i][j], S->L[i-1][j-1]+new_score, S->M[i-1][j-1]+new_score, S->U[i-1][j-1]+new_score);
			// LOW
			S->pointerL[i][j] = max3(&S->L[i][j], S->L[i-1][j]+EXTENSION, S->M[i-1][j]+GAP, -INFINITY);
			// UPP
			S->pointerU[i][j] = max3(&S->U[i][j], -INFINITY, S->M[i][j-1]+GAP, S->U[i][j-1]+EXTENSION);
		}
	}
	double max_score;
	int max_state = max3(&max_score, S->L[s1->l][s2->l], S->M[s1->l][s2->l], S->U[s1->l][s2->l]);
	trace_back(S, s1, s2, r1, r2, max_state);	
	destory_matrix(S);
	return max_score;
}

/* main function. */
int main(int argc, char *argv[]) {
	if((BLOSUM62 = load_BLOSUM62("test/BLOSUM62.txt")) == NULL) die("fail to load BLOSUM62 table at %s", "test/BLOSUM62.txt");
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
	printf("score=%f\n", align(ks1, ks2, r1, r2));
	printf("%s\n%s\n", r1->s, r2->s);
	kstring_destory(ks1);
	kstring_destory(ks2);
	kstring_destory(r1);
	kstring_destory(r2);
	return 0;
}