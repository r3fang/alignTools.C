all:
		$(CC) -g -O2 src/main.c src/utils.c src/kstring.c -o bin/alignTools -lz
		#$(CC) -g -O2 src/global_affine.c src/utils.c src/kstring.c -o bin/global_affine -lz
		#$(CC) -g -O2 src/local.c src/utils.c src/kstring.c -o bin/local -lz
		#$(CC) -g -O2 src/global.c src/utils.c src/kstring.c -o bin/global -lz
		#$(CC) -g -O2 src/local_affine.c src/utils.c src/kstring.c -o bin/local_affine -lz
		#$(CC) -g -O2 src/edit_dist.c src/utils.c src/kstring.c -o bin/edit_dist -lz
		#$(CC) -g -O2 src/fit.c src/utils.c src/kstring.c -o bin/fit -lz
		#$(CC) -g -O2 src/overlap.c src/utils.c src/kstring.c -o bin/overlap -lz
		#$(CC) -g -O2 src/fit_affine.c src/utils.c src/kstring.c -o bin/fit_affine -lz
		#$(CC) -g -O2 src/fit_affine_jump.c src/utils.c src/kstring.c -o bin/fit_affine_jump -lz		
clean:
		rm -f bin/*.dSYM
