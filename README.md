# LCP-Finder
Implement the skew algorithm to construct suffix array in O(n) time complexity and then find the longest common prefix array with Kauai's algorithm in linear time.

### Usage
First, Complile the main.cpp to executable file if needed.

Second, given a selected input DNA sequence as argument, the program will output the suffix array, LCP array and the execution time.
```
./main seq16.fasta

> Input sequence: ACGTGAGTAGCACTCG
> Suffix array:
 0 11  8  5 10 14  1 12 15  4  9  6  2  7 13  3
> LCP:
 2  1  2  0  1  2  1  0  1  1  1  2  0  1  1  0

Execution time: 21 ms
```

### Memo
The  ``n`` in  seq``n``.fasta file represent the number of nucleotide in the DNA sequence
