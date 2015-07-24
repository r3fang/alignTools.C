all:
		$(CC) -g -O2 src/main.c src/utils.c src/kstring.c -o bin/alignTools -lz
clean:
		rm -f bin/*.dSYM
