all:src/align.c src/utils.c src/utils.h src/kseq.h
		$(CC) -g -O2 src/align.c src/utils.c -o bin/align -lz

clean:
		rm -f bin/*.dSYM
		rm -r bin/*.dSYM
