#include <stdio.h>
#include <string.h>
#include "kstring.h"
#include "utils.h"
#include "overlap.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.7.23-r15"
#endif

int main_global_affine(int argc, char *argv[]);
int main_local_affine(int argc, char *argv[]);
int main_overlap(int argc, char *argv[]);
int main_fit_affine_jump(int argc, char *argv[]);
int main_edit_dist(int argc, char *argv[]);

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: alignTools (pairwise DNA sequence alignment)\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Rongxin Fang <r3fang@ucsd.edu>\n\n");
	fprintf(stderr, "Usage:   alignTools <command> [options]\n\n");
	fprintf(stderr, "Command: global     global (needle) alignment allows affine gap\n");
	fprintf(stderr, "         local      smith-waterman with affine gap\n");
	fprintf(stderr, "         fit        fit alingment allows affine gap plus jump state\n");
	fprintf(stderr, "         overlap    overlap alignment\n");
	fprintf(stderr, "         edit       edit distance\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	int i, ret;
	double t_real;
	kstring_t pg = {0,0,0};
	ksprintf(&pg, "@PG\tID:bwa\tPN:bwa\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);
	for (i = 1; i < argc; ++i) ksprintf(&pg, " %s", argv[i]);
	//bwa_pg = pg.s;
	if (argc < 2) return usage();
	else if (strcmp(argv[1], "global") == 0) ret = main_global_affine(argc-1, argv+1);
	else if (strcmp(argv[1], "local") == 0) ret = main_local_affine(argc-1, argv+1);
	else if (strcmp(argv[1], "fit") == 0) ret = main_fit_affine_jump(argc-1, argv+1);
	else if (strcmp(argv[1], "overlap") == 0) ret = main_overlap(argc-1, argv+1);
	else if (strcmp(argv[1], "edit") == 0) ret = main_edit_dist(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	//if (ret == 0) {
	//	fprintf(stderr, "[%s] Version: %s\n", __func__, PACKAGE_VERSION);
	//	fprintf(stderr, "[%s] CMD:", __func__);
	//	for (i = 0; i < argc; ++i)
	//		fprintf(stderr, " %s", argv[i]);
	//	fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());
	//}
	//free(bwa_pg);
	return ret;
}
