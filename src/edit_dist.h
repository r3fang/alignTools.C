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

void min3(int *res, int a1, int a2, int a3){
	*res = INFINITY;
	if(a1 < *res) *res = a1; 
	if(a2 < *res) *res = a2;
	if(a3 < *res) *res = a3;
}

int edit_dist(kstring_t *s1, kstring_t *s2){
	if(s1 == NULL || s2 == NULL) die("edit_dist: parameter error\n");
	size_t m   = s1->l + 1;
	size_t n   = s2->l + 1;
	matrix_t *S = create_matrix(m, n);
	for(i=0; i < S->m; i++) S->M[i][0] = i;
	for(j=0; j < S->n; j++) S->M[0][j] = j;
	
	size_t i, j, k, l;
	for(i = 1; i <= s1->l; i++){
		for(j = 1; j <= s2->l; j++){
	        int new_score = (strncmp(s1->s+(i-1), s2->s+(j-1), 1) == 0) ? 0 : 1;
			min3(&S->score[i][j], S->score[i][j-1] + 1, S->score[i-1][j-1] + new_score, S->score[i-1][j] + 1);
		}
	}
	destory_matrix(S);
	return S->score[s1->l][s2->l];
}

/* main function. */
int edit_dist_main(int argc, char *argv[]) {
	kstring_t *ks1, *ks2; 
	ks1 = mycalloc(1, kstring_t);
	ks2 = mycalloc(1, kstring_t);
	if (argc == 1) {
		fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
		return 1;
	}
	kstring_read(argv[1], ks1, ks2);
	if(ks1->s == NULL || ks2->s == NULL) die("fail to read sequence\n");
	printf("score=%d\n", edit_dist(ks1, ks2));
	kstring_destory(ks1);
	kstring_destory(ks2);
	return 0;
}