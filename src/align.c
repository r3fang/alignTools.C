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
#include "kstring.h"

KSEQ_INIT(gzFile, gzread);

#define GL_ERR_NONE 			0
#define GAP 					-1.0
#define MATCH 					2.0
#define MISMATCH 				-0.5

typedef enum {true, false} bool;

typedef struct {
  unsigned int m;
  unsigned int n;
  double **score;
} matrix_t;

matrix_t *create_matrix(size_t m, size_t n){
	size_t i, j; 
	matrix_t *S = mycalloc(1, matrix_t);
	S->m = m;
	S->n = n;
	S->score = mycalloc(m, double*);
	for (i = 0; i < m; i++) {
      S->score[i] = mycalloc(n, double);
    }
	return S;
}

int destory_matrix(matrix_t *S){
	if(S == NULL) die("destory_matrix: parameter error\n");
	int i, j;
	for(i = 0; i < S->m; i++){
		free(S->score[i]);
	}
	free(S);
	return GL_ERR_NONE;
}

double max3(double a1, double a2, double a3){
	double res = DBL_MIN;
    res = (a1 >= res) ? a1 : res;
    res = (a2 >= res) ? a2 : res;	
    res = (a3 >= res) ? a3 : res;	
	return res;
}

double global(kstring_t *s1, kstring_t *s2){
	size_t m   = s1->l + 1;
	size_t n   = s2->l + 1;
	matrix_t *S = create_matrix(m, n);
	size_t i, j, k, l;
	for(i=0; i < S->m; i++) S->score[i][0] = 0.0;
	for(j=0; j < S->n; j++) S->score[0][j] = 0.0;
	
	for(i = 1; i <= s1->l; i++){
		for(j = 1; j <= s2->l; j++){
	        double new_score = (strncmp(s1->s+(i-1), s2->s+(j-1), 1) == 0) ? MATCH : MISMATCH;
			S->score[i][j] = DBL_MIN;
			S->score[i][j] = max3(S->score[i-1][j-1] + new_score, S->score[i-1][j] + GAP, S->score[i][j-1] + GAP);
		}
	}
	double res = S->score[s1->l][s2->l];
	if(destory_matrix(S) != GL_ERR_NONE) die("smith_waterman: fail to destory matrix");
	return res;
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

void read_kstring(char* fname, kstring_t *str1, kstring_t *str2){
	gzFile fp;
	kseq_t *seq;
	fp = gzopen(fname, "r");
	seq = kseq_init(fp);
	int i, l;
	char **tmp = mycalloc(3, char*);
	i = 0;
	while((l=kseq_read(seq)) >= 0){
		if(i >= 2) die("input fasta file has more than 2 sequences");
		tmp[i++] = str_toupper(seq->seq.s);
	}
	if(tmp[0] == NULL || tmp[1] == NULL) die("read_kstring: fail to read sequence");
	(str1)->s = strdup(tmp[0]); (str1)->l = strlen((str1)->s);
	(str2)->s = strdup(tmp[1]); (str2)->l = strlen((str2)->s);
	for(; i >=0; i--){if(tmp[i]) free(tmp[i]);} free(tmp);
	kseq_destroy(seq);
	gzclose(fp);
}

void kstring_destory(kstring_t *ks){
	free(ks->s);
	free(ks);
}

/* main function. */
int main(int argc, char *argv[]) {	
	kstring_t *ks1, *ks2; 
	ks1 = mycalloc(1, kstring_t);
	ks2 = mycalloc(1, kstring_t);
	
	if (argc == 1) {
		fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
		return 1;
	}
	read_kstring(argv[1], ks1, ks2);
	if(ks1->s == NULL || ks2->s == NULL) die("fail to read sequence\n");
	printf("score=%f\n", global(ks1, ks2));
	kstring_destory(ks1);
	kstring_destory(ks2);
	return 0;
}