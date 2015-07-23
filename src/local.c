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
#define GAP                     -5.0

// POINTER STATE
#define LEFT                    100
#define DIAGONAL                200
#define RIGHT                   300
#define HOME                    400

typedef enum {true, false} bool;

/* DP matrix */
typedef struct {
  unsigned int m;
  unsigned int n;
  double **score;
  int  **pointer;
} matrix_t;

/* store scoring matrix */
scoring_matrix_t *BLOSUM62;

/*
 * create and initilize matrix
 */	
matrix_t *create_matrix(size_t m, size_t n){
	size_t i, j; 
	matrix_t *S = mycalloc(1, matrix_t);
	S->m = m;
	S->n = n;
	S->score = mycalloc(m, double*);
	for (i = 0; i < m; i++) {
      S->score[i] = mycalloc(n, double);
    }	
	
	S->pointer = mycalloc(m, int*);
	for (i = 0; i < m; i++) {
      S->pointer[i] = mycalloc(n, int);
    }
	// first row and first column initilized with 0's
	for(i=0; i < S->m; i++) S->score[i][0] = 0.0;
	for(j=0; j < S->n; j++) S->score[0][j] = 0.0;
	return S;
}

int destory_matrix(matrix_t *S){
	if(S == NULL) die("destory_matrix: parameter error\n");
	int i, j;
	for(i = 0; i < S->m; i++){
		free(S->score[i]);
	}
	for(i = 0; i < S->m; i++){
		free(S->pointer[i]);
	}
	free(S);
	return GL_ERR_NONE;
}

int max4(double *res, double a1, double a2, double a3, double a4){
	*res = -INFINITY;
	int state;
	if(a1 > *res){*res = a1; state = LEFT;}
	if(a2 > *res){*res = a2; state = DIAGONAL;}
	if(a3 > *res){*res = a3; state = RIGHT;}	
	if(a4 > *res){*res = a4; state = HOME;}
	return state;
}
/*
 * trace back starts from (i_max, j_max) and stops where M(i,j)=0.
 */
void trace_back(matrix_t *S, kstring_t *ks1, kstring_t *ks2, kstring_t *res_ks1, kstring_t *res_ks2, int i, int j){
	if(S == NULL || ks1 == NULL || ks2 == NULL || res_ks1 == NULL || res_ks2 == NULL) die("trace_back: parameter error");
	int m = 0; 
	while(S->score[i][j] != 0){
		switch(S->pointer[i][j]){
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
			case HOME: // go home baby
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
double align(kstring_t *s1, kstring_t *s2, kstring_t *r1, kstring_t *r2){
	if(s1 == NULL || s2 == NULL || r1 == NULL || r2 == NULL) die("global: parameter error\n");
	size_t m   = s1->l + 1;
	size_t n   = s2->l + 1;
	matrix_t *S = create_matrix(m, n);
	size_t i, j, k, l;
	double max_score = -INFINITY;
	int i_max, j_max;
	
	for(i = 1; i <= s1->l; i++){
		for(j = 1; j <= s2->l; j++){
			double new_score = match(s1->s[i-1], s2->s[j-1], BLOSUM62);
			S->pointer[i][j] = max4(&S->score[i][j], S->score[i][j-1] + GAP, S->score[i-1][j-1] + new_score, S->score[i-1][j] + GAP, 0.0);
			if(max_score < S->score[i][j]){
				max_score = S->score[i][j];
				i_max = i;
				j_max = j;
			}
		}
	}
	// find max value of S->score, can be anywhere in matrix
	// stop when we get a cell with 0
	trace_back(S, s1, s2, r1, r2, i_max, j_max);
	if(destory_matrix(S) != GL_ERR_NONE) die("smith_waterman: fail to destory matrix");
	return max_score;
}

/* main function. */
int main(int argc, char *argv[]) {
	if((BLOSUM62 = load_BLOSUM62("test/PAM250.txt")) == NULL) die("fail to load BLOSUM62 table at %s", "test/BLOSUM62.txt");
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
	printf("%f\n", align(ks1, ks2, r1, r2));
	printf("%s\n%s\n", r1->s, r2->s);
	kstring_destory(ks1);
	kstring_destory(ks2);
	kstring_destory(r1);
	kstring_destory(r2);
	return 0;
}