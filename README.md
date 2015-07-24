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
$ ./bin/alignTools 
Program: alignTools (pairwise DNA sequence alignment)
Version: 0.7.23-r15
Contact: Rongxin Fang <r3fang@ucsd.edu>

Usage:   alignTools <command> [options]

Command: gl         classic global alingment (needle)
         gla        global alignment with affine gap
         sw         classic smith-waterman alignment
         swa        smith-waterman with affine gap
         fit        fitting alignment
         fita       fitting alignment with affine gap
         fitaj      fitting alingment with affine gap plus jump state
         ov         overlap alignment
         ed         count edit distance
```

### Author
Rongxin Fang (r3fang@eng.ucsd.edu)