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
KSEQ_INIT(gzFile, gzread);


typedef struct {
  unsigned int m;
  double **score;
  char* bases;
} scoring_matrix_t;

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

static inline char 
*strrev(char *s){
	if(s == NULL) return NULL;
	int l = strlen(s);
	char *ss = strdup(s);
	free(s);
	s = mycalloc(l, char);
	int i; for(i=0; i<l; i++){
		s[i] = ss[l-i-1];
	}
	s[l] = '\0';
	return s;
}

/*
 * change string to upper case
 */
static inline char 
*str_toupper(char* s){
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
/*
 * destory kstring
 */
static inline void 
kstring_destory(kstring_t *ks){
	free(ks->s);
	free(ks);
}

static inline void 
kstring_read(char* fname, kstring_t *str1, kstring_t *str2){
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
