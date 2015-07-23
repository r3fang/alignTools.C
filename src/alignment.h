/*--------------------------------------------------------------------*/
/* alignment.h 		                                                  */
/* Author: Rongxin Fang                                               */
/* E-mail: r3fang@ucsd.edu                                            */
/* Date: 07-22-2015                                                   */
/* Basic functions/variables commonly used.                           */
/*--------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "utils.h"
#include "kseq.h"
#include "kstring.h"

typedef struct {
  unsigned int m;
  double **score;
  char* bases;
} scoring_matrix_t;

/* store scoring matrix */
scoring_matrix_t *BLOSUM62;

/* scoring */
static inline double 
match(char a, char b, scoring_matrix_t *S){
	if(S==NULL) die("scoring: input error");
	char *s_a, *s_b;
	int a_i, b_i;
	s_a = strchr (S->bases, a);
	s_b = strchr (S->bases, b);
	if(s_a == NULL || s_b == NULL) return -INFINITY;
	a_i = s_a - S->bases;
	b_i = s_b - S->bases;	
	return S->score[a_i][b_i];
}


/*
 * load the scoring matrix, can be BLOSUM62 or PAM250
 */
static inline scoring_matrix_t 
*load_BLOSUM62(char* fname){
	if(fname == NULL) die("load_BLOSUM62: input error");
	int i, j, n;
	int *fields;
	FILE *fp;
	scoring_matrix_t *S = mycalloc(1, scoring_matrix_t);
	kstring_t *buffer = mycalloc(1, kstring_t);
	buffer->s = mycalloc(4096, char);		
	if((fp=fopen(fname, "r")) == NULL) die("load_score_mat: %s not exists", fname);
	getline(&(buffer->s), &(buffer->l), fp);
	fields = ksplit(buffer, 0, &n);
	/* initilize the S*/
	S->m = n;
	S->score = mycalloc(n, double*);
	for(i=0; i<S->m; i++) S->score[i] = mycalloc(n, double);
	S->bases = mycalloc(n, char);
	for (j = 0; j < n; j++){S->bases[j] = *(buffer->s + fields[j]);}
	
	i=0;
	while ((getline(&(buffer->s), &(buffer->l), fp)) != -1) {
		fields = ksplit(buffer, 0, &n);
		for (j = 0; j < n; j++){
			S->score[i][j] = atof(buffer->s + fields[j]);			
		}
		i ++;
	}
	fclose(fp);
   	return S;
}