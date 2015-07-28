/*--------------------------------------------------------------------*/
/* alignment.h 		                                                  */
/* Author: Rongxin Fang                                               */
/* E-mail: r3fang@ucsd.edu                                            */
/* Date: 07-22-2015                                                   */
/* Basic functions/variables commonly used.                           */
/*--------------------------------------------------------------------*/
#ifndef _ALIGNMENT_
#define _ALIGNMENT_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <float.h>
#include <math.h>
#include <limits.h>		/* INT_MAX etc. */
#include <errno.h>
#include "zlib.h"
#include "kseq.h"
#include "kstring.h"

KSEQ_INIT(gzFile, gzread);
typedef enum { true, false } bool;

// POINTER STATE
#define LEFT                    100
#define DIAGONAL                200
#define RIGHT                   300
#define HOME                    400
#define LOW                     500
#define MID                     600
#define UPP                     700
#define JUMP                    800

// scoring matrix and pointer matrix
typedef struct {
  unsigned int m;
  unsigned int n;
  double **L;
  double **M;
  double **U;
  double **J;
  int  **pointerL;
  int  **pointerM;
  int  **pointerU;
  int  **pointerJ;
} matrix_t;

//for alignment allows jump state with junctions
typedef struct {
	size_t size;
	int *pos;
} junction_t;

//opt
typedef struct {
	int o; // gap open
	int e; // gap extension
	int m; // match
	int u; // unmatch
	int j; // jump penality
	bool s;
	junction_t sites;
} opt_t;


#define myalloc(n,type)	(type*)_myalloc((n)*sizeof(type))
#define mycalloc(n,type) (type*)_mycalloc(n,sizeof(type))

static inline void die (char *format, ...)
{
  va_list args ;

  va_start (args, format) ;
  fprintf (stderr, "FATAL ERROR: ") ;
  vfprintf (stderr, format, args) ;
  fprintf (stderr, "\n") ;
  va_end (args) ;
  exit (-1) ;
}

long int totalAllocated = 0 ;

static inline void *_mycalloc (long number, int size)
{
  void *p = (void*) calloc (number, size) ;
  if (!p) die ("mycalloc failure requesting %d of size %d bytes", number, size) ;
  totalAllocated += number*size ;
  return p ;
}

/* max of fix values */
static inline int 
max5(double *res, double a1, double a2, double a3, double a4, double a5){
	*res = -INFINITY;
	int state;
	if(a1 > *res){*res = a1; state = 0;}
	if(a2 > *res){*res = a2; state = 1;}
	if(a3 > *res){*res = a3; state = 2;}	
	if(a4 > *res){*res = a4; state = 3;}	
	if(a5 > *res){*res = a5; state = 4;}	
	return state;
}

static inline opt_t 
*init_opt(){
	opt_t *opt = mycalloc(1, opt_t);
	opt->o = -5.0;
	opt->e = -1.0;
	opt->m =  1.0;
	opt->u = -2.0;
	opt->j = -10.0;
	opt->s = false;
	opt->sites.size = 0;	
	opt->sites.pos = NULL;	
	return opt;
}

/*
 * create matrix, allocate memor
 */
static inline matrix_t 
*create_matrix(size_t m, size_t n){
	size_t i, j; 
	matrix_t *S = mycalloc(1, matrix_t);
	S->m = m;
	S->n = n;
	S->L = mycalloc(m, double*);
	S->M = mycalloc(m, double*);
	S->U = mycalloc(m, double*);
	S->J = mycalloc(m, double*);
	
	for (i = 0; i < m; i++) {
      S->M[i] = mycalloc(n, double);
      S->L[i] = mycalloc(n, double);
      S->U[i] = mycalloc(n, double);
      S->J[i] = mycalloc(n, double);
    }
	S->pointerM = mycalloc(m, int*);
	S->pointerU = mycalloc(m, int*);
	S->pointerL = mycalloc(m, int*);
	S->pointerJ = mycalloc(m, int*);
	for (i = 0; i < m; i++) {
       	S->pointerU[i] = mycalloc(n, int);
        S->pointerM[i] = mycalloc(n, int);
        S->pointerL[i] = mycalloc(n, int);
		S->pointerJ[i] = mycalloc(n, int);
    }
	return S;
}

/*
 * destory matrix
 */
static inline void 
destory_matrix(matrix_t *S){
	if(S == NULL) die("destory_matrix: parameter error\n");
	int i;
	for(i = 0; i < S->m; i++){
		if(S->L[i]) free(S->L[i]);
		if(S->M[i]) free(S->M[i]);
		if(S->U[i]) free(S->U[i]);
		if(S->J[i]) free(S->J[i]);
	}
	for(i = 0; i < S->m; i++){
		if(S->pointerL[i]) free(S->pointerL[i]);
		if(S->pointerM[i]) free(S->pointerM[i]);
		if(S->pointerU[i]) free(S->pointerU[i]);
		if(S->pointerJ[i]) free(S->pointerJ[i]);
	}
	free(S);
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

/*
 * read two sequences as str1 and str2 from the fasta file;
 * read junctions sites to opt->sites;
 * 
 */
static inline void 
kstring_read(char* fname, kstring_t *str1, kstring_t *str2, opt_t *opt){
	// input check
	if(fname == NULL || str1 == NULL || str2 == NULL || opt == NULL) 
		die("kstring_read: input error");
	
	// variables declarision
	int i, l; gzFile fp; kseq_t *seq;
	char **tmp_seq = mycalloc(3, char*);
	char **tmp_comment = mycalloc(3, char*);	
	// parser fasta
	fp = gzopen(fname, "r");
	seq = kseq_init(fp);
	if(fp == NULL || seq == NULL) die("Can't open %s\n", fname);
	
	i = 0; while((l=kseq_read(seq)) >= 0){
		if(i >= 2) die("input fasta file has more than 2 sequences");
		tmp_seq[i] = strdup(seq->seq.s);
		if(seq->comment.s) tmp_comment[i] = strdup(seq->comment.s);
		i++;
	}
	// read sequence
	if(tmp_seq[0] == NULL || tmp_seq[1] == NULL) die("read_kstring: fail to read sequence");
	(str1)->s = strdup(tmp_seq[0]); (str1)->l = strlen((str1)->s);
	(str2)->s = strdup(tmp_seq[1]); (str2)->l = strlen((str2)->s);
	// read the junctions sites if opt != NULL and opt->s==ture
	if(opt != NULL && opt->s == true){
		if(tmp_comment[1] == NULL) die("fail to read junction sites");
		kstring_t *tmp = mycalloc(1, kstring_t);
		tmp->s = strdup(tmp_comment[1]);
		tmp->l = strlen(tmp->s);
		int *fields, i, n;
		fields = ksplit(tmp, '|', &n);
		opt->sites.size = n;
		opt->sites.pos = mycalloc(n, int);
		for (i = 0; i < n; ++i) opt->sites.pos[i] = atoi(tmp->s + fields[i]);
		if(tmp) kstring_destory(tmp);
		if(fields) free(fields);
	}
	for(; i >=0; i--) if(tmp_seq[i]) free(tmp_seq[i]);	
	free(tmp_seq);
	if(tmp_comment) free(tmp_comment);
	if(seq) kseq_destroy(seq);
	gzclose(fp);
}

/*
 * check if an element in an array
 */
static inline bool 
isvalueinarray(int val, int *arr, int size){
    int i;
    for (i=0; i < size; i++) {
        if (arr[i] == val)
            return true;
    }
    return false;
}

/*
 * min value of three
 */
static inline void 
min3(double *res, double a1, double a2, double a3){
	*res = INFINITY;
	if(a1 < *res) *res = a1; 
	if(a2 < *res) *res = a2;
	if(a3 < *res) *res = a3;
}
/*--------------------------------------------------------------------*/
/* 
 * calculate edit distance
 */
static inline int 
edit_dist(kstring_t *s1, kstring_t *s2, opt_t *opt){
	if(s1 == NULL || s2 == NULL || opt == NULL) die("edit_dist: parameter error\n");
	double mismatch = opt->u;
	double match = 0;
	double gap = opt->o;
	size_t m   = s1->l + 1;
	size_t n   = s2->l + 1;
	matrix_t *S = create_matrix(m, n);
	size_t i, j, k, l;
	for(i=0; i < S->m; i++) S->M[i][0] = i;
	for(j=0; j < S->n; j++) S->M[0][j] = j;
	for(i = 1; i <= s1->l; i++){
		for(j = 1; j <= s2->l; j++){
			double new_score = ((s1->s[i-1] - s2->s[j-1]) == 0) ? match : mismatch;			
			min3(&S->M[i][j],
			      S->M[i][j-1] + 1, 
				  S->M[i-1][j-1] + new_score, 
				  S->M[i-1][j] + 1);
		}
	}
	int res = (int) S->M[s1->l][s2->l];
	destory_matrix(S);
	return res;
}

/* main function for edit dist */
static inline int 
main_edit_dist(int argc, char *argv[]) {	
	opt_t *opt = init_opt(); // initlize options with default settings
	int c;
	srand48(11);
	while ((c = getopt(argc, argv, "m:u:o:e")) >= 0) {
			switch (c) {
			case 'm': opt->m = atoi(optarg); break;
			case 'u': opt->u = atoi(optarg); break;
			case 'o': opt->o = atoi(optarg); break;
			case 'e': opt->e = atoi(optarg); break;
			default: return 1;
		}
	}
	if (optind + 1 > argc) {
		fprintf(stderr, "\n");
				fprintf(stderr, "Usage:   alignTools edit [options] <target.fa>\n\n");
				fprintf(stderr, "Options: -u INT   mismatch penalty [%d]\n", opt->u);
				fprintf(stderr, "         -o INT   gap penalty [%d]\n", opt->o);
				fprintf(stderr, "\n");
				return 1;
	}
	kstring_t *ks1, *ks2; 
	ks1 = mycalloc(1, kstring_t);
	ks2 = mycalloc(1, kstring_t);
	kstring_read(argv[argc-1], ks1, ks2, opt);
	if(ks1->s == NULL || ks2->s == NULL) die("fail to read sequence\n");
	printf("edit_distance=%d\n", edit_dist(ks1, ks2, opt));
	kstring_destory(ks1);
	kstring_destory(ks2);
	free(opt);
	return 0;
}
/*--------------------------------------------------------------------*/
/* 
 * global alignment allowing affine gap
 */
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
align_gla(kstring_t *s1, kstring_t *s2, kstring_t *r1, kstring_t *r2, opt_t *opt){
	if(s1 == NULL || s2 == NULL || r1 == NULL || r2 == NULL) die("align: parameter error\n");
	double mismatch = opt->u;
	double match = opt->m;
	double gap = opt->o;
	double extension = opt->e;
	
	size_t m   = s1->l + 1; size_t n   = s2->l + 1;
	matrix_t *S = create_matrix(m, n);
	// initlize DP matrix
	S->M[0][0] = 0.0;
	S->L[0][0] = S->U[0][0] = gap;
	// initlize 0 column
	int i, j;
	for(i=1; i<S->m; i++){
		S->L[i][0] = gap + extension*(i);
		S->M[i][0] = -INFINITY;
		S->U[i][0] = -INFINITY;
	}
	for(j=1; j<S->n; j++){
		S->L[0][j] = -INFINITY;
		S->M[0][j] = -INFINITY;
		S->U[0][j] = gap + extension*(j);
	}
	//-------------------------------
	double new_score;
	int idx;
	// recurrance relation
	for(i=1; i<=s1->l; i++){
		for(j=1; j<=s2->l; j++){
			// MID
			new_score = ((s1->s[i-1] - s2->s[j-1]) == 0) ? match : mismatch;
			//new_score = match(s1->s[i-1], s2->s[j-1], BLOSUM62);
			idx = max5(&S->M[i][j], S->L[i-1][j-1]+new_score, S->M[i-1][j-1]+new_score, S->U[i-1][j-1]+new_score, -INFINITY, -INFINITY);
			if(idx==0) S->pointerM[i][j] = LOW;
			if(idx==1) S->pointerM[i][j] = MID;
			if(idx==2) S->pointerM[i][j] = UPP;
			// LOW
			idx = max5(&S->L[i][j], S->L[i-1][j]+extension, S->M[i-1][j]+gap, -INFINITY, -INFINITY, -INFINITY);
			if(idx==0) S->pointerL[i][j] = LOW;
			if(idx==1) S->pointerL[i][j] = MID;			
			// UPP
			idx = max5(&S->U[i][j], -INFINITY, S->M[i][j-1]+gap, S->U[i][j-1]+extension, -INFINITY, -INFINITY);
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

/* main function for global alignment. */
static inline int 
main_global_affine(int argc, char *argv[]) {
	opt_t *opt = init_opt(); // initlize options with default settings
	int c;
	srand48(11);
	while ((c = getopt(argc, argv, "m:u:o:e:j:s")) >= 0) {
			switch (c) {
			case 'm': opt->m = atoi(optarg); break;
			case 'u': opt->u = atoi(optarg); break;
			case 'o': opt->o = atoi(optarg); break;
			case 'e': opt->e = atoi(optarg); break;
			default: return 1;
		}
	}
	if (optind + 1 > argc) {
		fprintf(stderr, "\n");
				fprintf(stderr, "Usage:   alignTools global [options] <target.fa>\n\n");
				fprintf(stderr, "Options: -m INT   score for a match [%d]\n", opt->m);
				fprintf(stderr, "         -u INT   mismatch penalty [%d]\n", opt->u);
				fprintf(stderr, "         -o INT   gap open penalty [%d]\n", opt->o);
				fprintf(stderr, "         -e INT   gap extension penalty [%d]\n", opt->e);
				fprintf(stderr, "\n");
				return 1;
	}
	kstring_t *ks1, *ks2; 
	ks1 = mycalloc(1, kstring_t);
	ks2 = mycalloc(1, kstring_t);
	kstring_read(argv[argc-1], ks1, ks2, opt);
	if(ks1->s == NULL || ks2->s == NULL) die("fail to read sequence\n");
	kstring_t *r1 = mycalloc(1, kstring_t);
	kstring_t *r2 = mycalloc(1, kstring_t);
	r1->s = mycalloc(ks1->l + ks2->l, char);
	r2->s = mycalloc(ks1->l + ks2->l, char);
	printf("score=%f\n", align_gla(ks1, ks2, r1, r2, opt));
	printf("%s\n%s\n", r1->s, r2->s);
	free(opt);
	kstring_destory(ks1);
	kstring_destory(ks2);
	kstring_destory(r1);
	kstring_destory(r2);
	return 0;
}

/*--------------------------------------------------------------------*/
/* fit alignment for fit alignment*/
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
align_fit_affine_jump(kstring_t *s1, kstring_t *s2, kstring_t *r1, kstring_t *r2, opt_t *opt){
	if(s1 == NULL || s2 == NULL || r1 == NULL || r2 == NULL || opt == NULL) die("align: parameter error\n");
	if(s1->l > s2->l) die("first sequence must be shorter than the second to do fitting alignment"); 
	size_t m   = s1->l + 1; size_t n   = s2->l + 1;
	matrix_t *S = create_matrix(m, n);
	// copy alignment parameter
	junction_t junctions = opt->sites;
	double match = opt->m;
	double mismatch = opt->u;
	double gap = opt->o;
	double extension = opt->e;
	double jump_penality = opt->j;
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
			new_score = ((s1->s[i-1] - s2->s[j-1]) == 0) ? match : mismatch;
			//new_score = (strnicmp(s1->s+(i-1), s2->s+(j-1), 1) == 0) ? match : mismatch;
			if(opt->s == true){
				idx = max5(&S->M[i][j], S->L[i-1][j-1]+new_score, S->M[i-1][j-1]+new_score, S->U[i-1][j-1]+new_score, S->J[i-1][j-1]+new_score, -INFINITY);
				if(idx == 0) S->pointerM[i][j]=LOW;
				if(idx == 1) S->pointerM[i][j]=MID;
				if(idx == 2) S->pointerM[i][j]=UPP;
				if(idx == 3) S->pointerM[i][j]=JUMP;			 				
			}else{
				idx = max5(&S->M[i][j], S->L[i-1][j-1]+new_score, S->M[i-1][j-1]+new_score, S->U[i-1][j-1]+new_score, -INFINITY, -INFINITY);
				if(idx == 0) S->pointerM[i][j]=LOW;
				if(idx == 1) S->pointerM[i][j]=MID;
				if(idx == 2) S->pointerM[i][j]=UPP;				
			}
			
			// LOW
			idx = max5(&S->L[i][j], S->L[i-1][j]+extension, S->M[i-1][j]+gap, -INFINITY, -INFINITY, -INFINITY);
			if(idx == 0) S->pointerL[i][j]=LOW;
			if(idx == 1) S->pointerL[i][j]=MID;
			
			// UPP
			idx = max5(&S->U[i][j], -INFINITY, S->M[i][j-1]+gap, S->U[i][j-1]+extension, -INFINITY, -INFINITY);
			if(idx == 1) S->pointerU[i][j]=MID;
			if(idx == 2) S->pointerU[i][j]=UPP;
			
			// JUMP only allowed going to JUMP state at junction sites
			if(opt->s == true){
				if(isvalueinarray(j-1, junctions.pos, junctions.size)){
					idx = max5(&S->J[i][j], -INFINITY, S->M[i][j-1]+jump_penality, -INFINITY, S->J[i][j-1], -INFINITY);
					if(idx == 1) S->pointerJ[i][j] = MID;			
					if(idx == 3) S->pointerJ[i][j] = JUMP;			
				}else{
					idx = max5(&S->J[i][j], -INFINITY, -INFINITY, -INFINITY, S->J[i][j-1], -INFINITY);				
					if(idx == 3) S->pointerJ[i][j] = JUMP;
				}	
			}
		}
	}
	
	// find trace-back start point
	// NOTE: ALWAYS STARTS TRACING BACK FROM MID OR LOW
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
	opt_t *opt = init_opt(); // initlize options with default settings
	int c;
	srand48(11);
	while ((c = getopt(argc, argv, "m:u:o:e:j:s")) >= 0) {
			switch (c) {
			case 'm': opt->m = atoi(optarg); break;
			case 'u': opt->u = atoi(optarg); break;
			case 'o': opt->o = atoi(optarg); break;
			case 'e': opt->e = atoi(optarg); break;
			case 'j': opt->j = atoi(optarg); break;
			case 's': opt->s = true; break;
			default: return 1;
		}
	}
	if (optind + 1 > argc) {
		fprintf(stderr, "\n");
				fprintf(stderr, "Usage:   alignTools fit [options] <target.fa>\n\n");
				fprintf(stderr, "Options: -m INT   score for a match [%d]\n", opt->m);
				fprintf(stderr, "         -u INT   mismatch penalty [%d]\n", opt->u);
				fprintf(stderr, "         -o INT   gap open penalty [%d]\n", opt->o);
				fprintf(stderr, "         -e INT   gap extension penalty [%d]\n", opt->e);
				fprintf(stderr, "         -j INT   jump penality [%d]\n", opt->j);
				fprintf(stderr, "         -s       weather jump state include\n");
				fprintf(stderr, "\n");
				return 1;
	}
	kstring_t *ks1, *ks2; 
	ks1 = mycalloc(1, kstring_t);
	ks2 = mycalloc(1, kstring_t);
	kstring_read(argv[argc-1], ks1, ks2, opt);
	if(ks1->s == NULL || ks2->s == NULL) die("fail to read sequence\n");
	if(ks1->l > ks2->l) die("first sequence must be shorter than the second\n");
	kstring_t *r1 = mycalloc(1, kstring_t);
	kstring_t *r2 = mycalloc(1, kstring_t);
	r1->s = mycalloc(ks1->l + ks2->l, char);
	r2->s = mycalloc(ks1->l + ks2->l, char);
	printf("score=%f\n", align_fit_affine_jump(ks1, ks2, r1, r2, opt));
	printf("%s\n%s\n", r1->s, r2->s);
	kstring_destory(ks1);
	kstring_destory(ks2);
	kstring_destory(r1);
	kstring_destory(r2);
	free(opt);
	return 0;
}


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
align_local_affine(kstring_t *s1, kstring_t *s2, kstring_t *r1, kstring_t *r2, opt_t *opt){
	if(s1 == NULL || s2 == NULL || r1 == NULL || r2 == NULL) die("align: parameter error\n");
	double mismatch = opt->u;
	double match = opt->m;
	double gap = opt->o;
	double extension = opt->e;
	
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
			new_score = ((s1->s[i-1] - s2->s[j-1]) == 0) ? match : mismatch;
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
			idx = max5(&S->L[i][j], S->L[i-1][j]+extension, S->M[i-1][j]+gap, -INFINITY, -INFINITY, -INFINITY);
			if(idx==0) S->pointerL[i][j] = LOW;
			if(idx==1) S->pointerL[i][j] = MID;
			// UPP
			idx = max5(&S->U[i][j], -INFINITY, S->M[i][j-1]+gap, S->U[i][j-1]+extension, -INFINITY, -INFINITY);
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
	opt_t *opt = init_opt(); // initlize options with default settings
	int c;
	srand48(11);
	while ((c = getopt(argc, argv, "m:u:o:e:j:s")) >= 0) {
			switch (c) {
			case 'm': opt->m = atoi(optarg); break;
			case 'u': opt->u = atoi(optarg); break;
			case 'o': opt->o = atoi(optarg); break;
			case 'e': opt->e = atoi(optarg); break;
			default: return 1;
		}
	}
	if (optind + 1 > argc) {
		fprintf(stderr, "\n");
				fprintf(stderr, "Usage:   alignTools local [options] <target.fa>\n\n");
				fprintf(stderr, "Options: -m INT   score for a match [%d]\n", opt->m);
				fprintf(stderr, "         -u INT   mismatch penalty [%d]\n", opt->u);
				fprintf(stderr, "         -o INT   gap open penalty [%d]\n", opt->o);
				fprintf(stderr, "         -e INT   gap extension penalty [%d]\n", opt->e);
				fprintf(stderr, "\n");
				return 1;
	}
	
	kstring_t *ks1, *ks2; 
	ks1 = mycalloc(1, kstring_t);
	ks2 = mycalloc(1, kstring_t);
	kstring_read(argv[argc-1], ks1, ks2, opt);
	if(ks1->s == NULL || ks2->s == NULL) die("fail to read sequence\n");
	kstring_t *r1 = mycalloc(1, kstring_t);
	kstring_t *r2 = mycalloc(1, kstring_t);
	r1->s = mycalloc(ks1->l + ks2->l, char);
	r2->s = mycalloc(ks1->l + ks2->l, char);
	printf("score=%f\n", align_local_affine(ks1, ks2, r1, r2, opt));
	printf("%s\n%s\n", r1->s, r2->s);
	kstring_destory(ks1);
	kstring_destory(ks2);
	kstring_destory(r1);
	kstring_destory(r2);
	return 0;
}
/*
 * trace back starts from (i_max, j_max) and stops where M(i,j)=0.
 */
static inline void trace_back_overlap(matrix_t *S, kstring_t *ks1, kstring_t *ks2, kstring_t *res_ks1, kstring_t *res_ks2, int i, int j){
	if(S == NULL || ks1 == NULL || ks2 == NULL || res_ks1 == NULL || res_ks2 == NULL) die("trace_back: parameter error");
	int m = 0; 
	while(j>0){
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
static inline double align_overlap(kstring_t *s1, kstring_t *s2, kstring_t *r1, kstring_t *r2, opt_t *opt){
	if(s1 == NULL || s2 == NULL || r1 == NULL || r2 == NULL) die("align_overlap: parameter error\n");
	double mismatch = opt->u;
	double match = opt->m;
	double gap = opt->o;
	double extension = opt->e;
	size_t m   = s1->l + 1;
	size_t n   = s2->l + 1;
	size_t i, j, k, l;
	matrix_t *S = create_matrix(m, n);
	// first row and first column initilized with 0's
	for(j=0; j < S->n; j++) S->M[0][j] = -INFINITY;
	for(i=0; i < S->m; i++) S->M[i][0] = 0.0;
	int idx;
	for(i = 1; i <= s1->l; i++){
		for(j = 1; j <= s2->l; j++){
			//double new_score = match(s1->s[i-1], s2->s[j-1], BLOSUM62);
			double new_score = ((s1->s[i-1] - s2->s[j-1]) == 0) ? match : mismatch;			
			idx = max5(&S->M[i][j], S->M[i][j-1] + gap, S->M[i-1][j-1] + new_score, S->M[i-1][j] + gap, -INFINITY, -INFINITY);
			if(idx==0) S->pointerM[i][j] = LEFT;
			if(idx==1) S->pointerM[i][j] = DIAGONAL;
			if(idx==2) S->pointerM[i][j] = RIGHT;			
		}
	}
	// find max value of on the bottom column of S->score and starts tracing back from there
	double max_score = -INFINITY;
	int i_max, j_max;
	i_max = s1->l;
	for(j=0; j<s2->l; j++){
		if(max_score < S->M[i_max][j]){
			max_score = S->M[i_max][j];
			j_max = j;
		}
	}
	// stop when we get a cell with 0
	trace_back_overlap(S, s1, s2, r1, r2, i_max, j_max);
	destory_matrix(S);
	return max_score;
}

/* main function for overlap alignment. */
static inline int main_overlap(int argc, char *argv[]) {
	opt_t *opt = init_opt(); // initlize options with default settings
	int c;
	srand48(11);
	while ((c = getopt(argc, argv, "m:u:o:e:j:s")) >= 0) {
			switch (c) {
			case 'm': opt->m = atoi(optarg); break;
			case 'u': opt->u = atoi(optarg); break;
			case 'o': opt->o = atoi(optarg); break;
			case 'e': opt->e = atoi(optarg); break;
			default: return 1;
		}
	}
	if (optind + 1 > argc) {
		fprintf(stderr, "\n");
				fprintf(stderr, "Usage:   alignTools overlap [options] <target.fa>\n\n");
				fprintf(stderr, "Options: -m INT   score for a match [%d]\n", opt->m);
				fprintf(stderr, "         -u INT   mismatch penalty [%d]\n", opt->u);
				fprintf(stderr, "         -o INT   gap open penalty [%d]\n", opt->o);
				fprintf(stderr, "         -e INT   gap extension penalty [%d]\n", opt->e);
				fprintf(stderr, "\n");
				return 1;
	}
	
	kstring_t *ks1, *ks2;
	ks1 = mycalloc(1, kstring_t);
	ks2 = mycalloc(1, kstring_t);
	kstring_read(argv[1], ks1, ks2, opt);
	if(ks1->s == NULL || ks2->s == NULL) die("fail to read sequence\n");
	kstring_t *r1 = mycalloc(1, kstring_t);
	kstring_t *r2 = mycalloc(1, kstring_t);
	r1->s = mycalloc(ks1->l + ks2->l, char);
	r2->s = mycalloc(ks1->l + ks2->l, char);
	printf("%f\n", align_overlap(ks1, ks2, r1, r2, opt));
	printf("%s\n%s\n", r1->s, r2->s);
	free(opt);
	kstring_destory(ks1);
	kstring_destory(ks2);
	kstring_destory(r1);
	kstring_destory(r2);
	return 0;
}
#endif