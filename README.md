# alignTools

alignTools is a collection of C-implemented tools for pairwise DNA sequence alignment. It includes following alignment features:

  - global alignment
  - local alingment
  - fit alignment (jump state)
  - overlap alignment 
  - edit distance

More features will be added in the near future.

## Version
0.7.23-r15

## Get Started

  - Install

```
$ git clone git@github.com:r3fang/alignTools.git
$ cd alignTools
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

  - global alingment

```
$./bin/alignTools global
Usage:   alignTools global [options] <target.fa>

Options: -m INT   score for a match [1]
         -u INT   mismatch penalty [-2]
         -o INT   gap open penalty [-5]
         -e INT   gap extension penalty [-1]

$./bin/alignTools global -m 1 -u -1 -o -4 -e -1 test/test_global.fa
```

  - local alingment

```
$./bin/alignTools local
Usage:   alignTools local [options] <target.fa>

Options: -m INT   score for a match [1]
         -u INT   mismatch penalty [-2]
         -o INT   gap open penalty [-5]
         -e INT   gap extension penalty [-1]

$./bin/alignTools local -m 2 -u -2 -o -5 -e -2 test/test_local.fa
```

  - fit alingment

```
$ ./bin/alignTools fit 
Usage:   alignTools fit [options] <target.fa>

Options: -m INT   score for a match [1]
         -u INT   mismatch penalty [-2]
         -o INT   gap open penalty [-5]
         -e INT   gap extension penalty [-1]
         -j INT   jump penality [-10]
         -s       weather jump state included

$./bin/alignTools fit -m 2 -u -2 -s test/test_fit.fa
```

  - overlap alignment

```
$./bin/alignTools overlap 
Usage:   alignTools overlap [options] <target.fa>

Options: -m INT   score for a match [1]
         -u INT   mismatch penalty [-2]
         -o INT   gap open penalty [-5]
         -e INT   gap extension penalty [-1]

$./bin/alignTools overlap test/test_overlap.fa
```

  - edit distance

```
$./bin/alignTools edit
Usage:   alignTools edit [options] <target.fa>

Options: -u INT   mismatch penalty [-2]
         -o INT   gap penalty [-5]

$./bin/alignTools edit -u 1 -o 2 test/test_edit.fa
```

## Author
Rongxin Fang (r3fang@eng.ucsd.edu)
