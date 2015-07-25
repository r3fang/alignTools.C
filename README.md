# AlignTools

AlignTools is a collection of C-implemented pairwise DNA sequence alignment tools. It includes following alignment methods:

  - needle (global alignment) 
  - smith-waterman (local alingment)
  - global/local alignment allows affine gap
  - edit distance
  - overlap alignment
  - fitting alignment allows affine gap plus jump state

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

Command: gl         classic global alingment (needle)
         gla        global alignment with affine gap
         sw         classic smith-waterman alignment
         swa        smith-waterman with affine gap
         fit        fitting alingment allows affine gap plus jump state
         ov         overlap alignment
         ed         count edit distance
```
fitting alignment example:
```
$ ./bin/alignTools fit 
Usage:   alignTools fit [options] <target.fa>

Options: -m INT   score for a match [1]
         -u INT   mismatch penalty [-2]
         -o INT   gap open penalty [-5]
         -e INT   gap extension penalty [-1]
         -j INT   jump penality [-10]
         -s       weather jump state included

$./bin/alignTools fit -m 2 -u -2 -s test/test_fit_affine_jump.fa
score=172.000000
AAACGCCAGACATA----GAACTAAATTCGCCCCCGAGCAACC------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAGACTASFSADFTCAGTCCGGAAAATTTACTGTTTTATGCGGCACTGTTGCCGCTCA
AAACGCCAGACATAGGTCGAACTAAATTCGCCCCCGAGCAACCGTGCGTGGAGAGTGAAGTTGGGTTATATAGTTCGGACTGTGTTTGGTGATTCGGTACTAGGTTGCGGTGTCACTTGATAACTCATCCCTACATTTTTGTTTTTAGAGGCTCTCGCCGAAGCGATTCATCACCAGTTACGGAAGCACTGTGGAGCGAGCACAACAAGGGTCGGGAGCGGCTTTAGCCACTACGAACGCATGAGGGCGATTTCGACAGCGCAGTGTATCGCGGACCACAGACGTATGAAATACCCACAGCAGACATGCGACATTCGCGAGGTTGGCAAAGTCACTTATTCCAGTGGGTGCATCCCTTCAGGACCCACTCAGATCTCGCAGGCCACTCATCCACTATACTGCCACTAAGTCTCGGTACATAAGGGGCATATCATCATATACCTACGCGTGCCATTAAATTATGACTACCTGGGATCGAGGACCTTCGTGTGAAACAGACTAGGCACAACGCCAGACTASFSADFTCAGTCCGGAAAATTTACTGTTTTATGCGGCACTGTTGCCGCTCA
```
### Author
Rongxin Fang (r3fang@eng.ucsd.edu)