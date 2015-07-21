/*--------------------------------------------------------------------*/
/* GlobalJumpAligner.c                                                */
/* Author: Rongxin Fang                                               */
/* E-mail: r3fang@ucsd.edu                                            */
/* Predict Gene Fusion by given fastq files.                          */
/*--------------------------------------------------------------------*/
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread);


#define GL_ERR_NONE 			0
#define GAP 					-1.0
#define MATCH 					2.0
#define MISMATCH 				-0.5

typedef enum {true, false} bool;

typedef struct {
  char *seq1;
  unsigned int len1;
  char *seq2;
  unsigned int len2;
} seq_pair;
typedef seq_pair *seq_pair_t;

typedef struct {
  unsigned int m;
  unsigned int n;
  double **mat;
} matrix;
typedef matrix *matrix_t;

matrix_t create_matrix(size_t m, size_t n){
	size_t i, j; 
	matrix_t S = mycalloc(1, matrix);
	S->m = m;
	S->n = n;
	S->mat = mycalloc(m, double*);
	for (i = 0; i < m; i++) {
      S->mat[i] = mycalloc(n, double);
    }
	return S;
}

int destory_matrix(matrix_t S){
	if(S == NULL) die("destory_matrix: parameter error\n");
	int i, j;
	for(i = 0; i < S->m; i++){
		free(S->mat[i]);
	}
	free(S);
	return GL_ERR_NONE;
}

int destory_seq_pair(seq_pair_t p){
	if(p == NULL) die("destory_seq_pair: parameter error\n");
	free(p->seq1);
	free(p->seq2);
	free(p);
	return GL_ERR_NONE;
}

double max3(double a1, double a2, double a3){
	double res = DBL_MIN;
    res = (a1 >= res) ? a1 : res;
    res = (a2 >= res) ? a2 : res;	
    res = (a3 >= res) ? a3 : res;	
	return res;
}

double global(seq_pair *problem, bool local){
	size_t m   = problem->len1 + 1;
	size_t n   = problem->len2 + 1;
	matrix_t S = create_matrix(m, n);
	size_t i, j, k, l;
	for(i=0; i < S->m; i++) S->mat[i][0] = 0.0;
	for(j=0; j < S->n; j++) S->mat[0][j] = 0.0;
	
	for(i = 1; i <= problem->len1; i++){
		for(j = 1; j <= problem->len2; j++){
	        double new_score = (strncmp(problem->seq1+(i-1), problem->seq2+(j-1), 1) == 0) ? MATCH : MISMATCH;
			S->mat[i][j] = DBL_MIN;
			S->mat[i][j] = max3(S->mat[i-1][j-1] + new_score, S->mat[i-1][j] + GAP, S->mat[i][j-1] + GAP);
		}
	}
	double res = S->mat[problem->len1][problem->len2];
	if(destory_matrix(S) != GL_ERR_NONE) die("smith_waterman: fail to destory matrix");
	return res;
}

seq_pair_t init_pair_t(char* str1, char* str2){
	if(str1 == NULL || str2 == NULL) die("init_pair_t: parameter error\n");
    seq_pair_t problem = mycalloc(1, seq_pair);
	problem->seq1 = mycalloc(strlen(str1), char);
	problem->seq2 = mycalloc(strlen(str2), char);
	problem->seq1 = strdup(str1);
	problem->len1 = strlen(problem->seq1);
	problem->seq2 = strdup(str2);
	problem->len2 = strlen(problem->seq2);
	return problem;
}

char* str_toupper(char* s){
	char *r = mycalloc(strlen(s), char);
	int i = 0;
	char c;
	while(s[i])
	{
		r[i] = toupper(s[i]);
		i++;
	}
	r[strlen(s)] = '\0';
	return r;
}

/* main function. */
int main(int argc, char *argv[]) {	
	gzFile fp;
	kseq_t *seq;
	char *str1, *str2; 
	int l;
	if (argc == 1) {
		fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
		return 1;
	}
	fp = gzopen(argv[1], "r");
	seq = kseq_init(fp);
	if((l=kseq_read(seq)) >= 0){
		str1 = str_toupper(seq->seq.s);
	}
	if((l=kseq_read(seq) >= 0)){
		str2 = str_toupper(seq->seq.s);
	}
	if(str1 == NULL || str2 == NULL) die("fail to read sequence\n");
    printf("%s\t%s\n", str1, str2);
	seq_pair_t problem = init_pair_t(str1, str2);
	printf("score=%f\n", global(problem, false));
	if(destory_seq_pair(problem) != GL_ERR_NONE) die("destory_seq_pair fails\n");
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}