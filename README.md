# AlignTools

AlignTools is a collection of C-implemented pairwise DNA sequence alignment tools. It includes following alignment methods:

  - needle (global alignment) 
  - smith-waterman (local alingment)
  - global/local alignment with affine gap
  - fitting alignment
  - fitting alingment with affine gap
  - fitting alignment with affine gap and jump state

More features may be added in the near future.
### Version
1.0.1

### Get Started

Install from github:

```sh
$ git clone git@github.com:r3fang/alignTools.git
$ cd AlignTools
$ make
$ ./bin/global test/test.fa
$ ./bin/fit_affine_jump test/test_fit_affine_jump.fa 
```

### Author
Rongxin Fang (r3fang@eng.ucsd.edu)