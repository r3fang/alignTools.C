# AlignTools

AlignTools is a collection of C-implemented pairwise DNA sequence alignment tools. It includes following alignment methods:

  - needle (global alignment) allows affine gap
  - smith-waterman (local alingment) allows affine gap
  - fit alignment allows affine gap plus jump state 
    (often used for rna-seq alignment with splicing junctions)
  - overlap alignment
  - edit distance

More features will be added in the near future.
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

Command: global     global (needle) alignment allows affine gap
         local      smith-waterman alignment with affine gap
         fit        fit alingment allows affine gap plus jump state
         overlap    overlap alignment
         edit       edit distance
```

  - fit alingment allows affine gap plus jump state

```
$ ./bin/alignTools fit 
Usage:   alignTools fit [options] <target.fa>

Options: -m INT   score for a match [1]
         -u INT   mismatch penalty [-2]
         -o INT   gap open penalty [-5]
         -e INT   gap extension penalty [-1]
         -j INT   jump penality [-10]
         -s       weather jump state included

$./bin/alignTools fit -m 2 -u -2 -s test/ABI1.fa
```

  - global alingment with affine gap

```
$./bin/alignTools global test/test_global_affine.fa

```

  - local alingment with affine gap

```
$./bin/alignTools local test/test_local_affine.fa

```

  - edit distance

```
$./bin/alignTools edit test/test.fa

```

### Author
Rongxin Fang (r3fang@eng.ucsd.edu)