all:src/global_affine_gap.c src/utils.c src/utils.h src/kseq.h src/kstring.c src/kstring.h
		$(CC) -g -O2 src/global_affine_gap.c src/utils.c src/kstring.c -o bin/align -lz

clean:
		rm -f bin/*.dSYM
