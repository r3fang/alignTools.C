/*--------------------------------------------------------------------*/
/* edit_dist.h 		                                                  */
/* Author: Rongxin Fang                                               */
/* E-mail: r3fang@ucsd.edu                                            */
/* Date: 07-22-2015                                                   */
/* Pair-wise edit distance.                                           */
/* initilize S(i,j):                                                  */
/* S(i, 0) = i; S(0, j) = j                                           */
/* reccurrance relations:                                             */
/* S(i,j) = min{S(i-1, j-1)+1 (if x!=y), S(i-1, j)+1, S(i, j-1)+1}    */
/* no trace back needed.                                              */
/*--------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "utils.h"     // die and mycalloc
#include "kseq.h"      // fasta parser
#include "kstring.h"   // kstring_t type
#include "alignment.h"

static inline void 
min3(double *res, double a1, double a2, double a3){
	*res = INFINITY;
	if(a1 < *res) *res = a1; 
	if(a2 < *res) *res = a2;
	if(a3 < *res) *res = a3;
}
/* calculate edit distance*/
static inline int 
edit_dist(kstring_t *s1, kstring_t *s2){
	if(s1 == NULL || s2 == NULL) die("edit_dist: parameter error\n");
	size_t m   = s1->l + 1;
	size_t n   = s2->l + 1;
	matrix_t *S = create_matrix(m, n);
	size_t i, j, k, l;
	for(i=0; i < S->m; i++) S->M[i][0] = i;
	for(j=0; j < S->n; j++) S->M[0][j] = j;
	for(i = 1; i <= s1->l; i++){
		for(j = 1; j <= s2->l; j++){
			int new_score = ((s1->s[i-1] - s2->s[j-1]) == 0) ? 0 : 1;			
			min3(&S->M[i][j], S->M[i][j-1] + 1, S->M[i-1][j-1] + new_score, S->M[i-1][j] + 1);
		}
	}
	int res = (int) S->M[s1->l][s2->l];
	destory_matrix(S);
	return res;
}

/* main function. */
static inline int 
main_edit_dist(int argc, char *argv[]) {
	kstring_t *ks1, *ks2; 
	opt_t *opt = mycalloc(1, opt_t);
	ks1 = mycalloc(1, kstring_t);
	ks2 = mycalloc(1, kstring_t);
	if (argc == 1) {
		fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
		return 1;
	}
	kstring_read(argv[1], ks1, ks2, opt);
	if(ks1->s == NULL || ks2->s == NULL) die("fail to read sequence\n");
	printf("edit_distance=%d\n", edit_dist(ks1, ks2));
	kstring_destory(ks1);
	kstring_destory(ks2);
	free(opt);
	return 0;
}