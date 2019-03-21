#!/usr/bin/env python

##############################################################################
#
# Name: tran.py
# Date: 2019-03-21
# Author: kimoton
# Description:
# Problem
# For DNA strings s1 and s2 having the same length, their
# transition/transversion ratio R(s1,s2) is the ratio of the total number of
# transitions to the total number of transversions, where symbol substitutions
# are inferred from mismatched corresponding symbols as when calculating
# Hamming distance (see "Counting Point Mutation").
#
# Given: Two DNA strings s1 and s2 of equal length (at most 1 kbp).

# Return: The transition/transversion ratio R(s1,s2).

# Sample Dataset
# >Rosalind_0209
# GCAACGCACAACGAAAACCCTTAGGGACTGGATTATTTCGTGATCGTTGTAGTTATTGGA
# AGTACGGGCATCAACCCAGTT
# >Rosalind_2200
# TTATCTGACAAAGAAAGCCGTCAACGGCTGGATAATTTCGCGATCGTGCTGGTTACTGGC
# GGTACGAGTGTTCCTTTGGGT

# Sample Output
# 1.21428571429
##############################################################################

import sys
from Bio import SeqIO

TRANSITION = [("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")]
TRANSVERSION = [("A", "C"), ("C", "A"), ("C", "G"), ("G", "C"),
                ("A", "T"), ("T", "A"), ("T", "G"), ("G", "T")]


def get_seq(fasta):
    records = []
    with open(fasta, "r") as rf:
        for record in SeqIO.parse(rf, "fasta"):
            records.append(str(record.seq))
    return records


def count_transversion(seq1, seq2):
    count = 0
    for i, j in zip(seq1, seq2):
        if (i, j) in TRANSVERSION:
            count += 1
    return count


def count_transition(seq1, seq2):
    count = 0
    for i, j in zip(seq1, seq2):
        if (i, j) in TRANSITION:
            count += 1
    return count


def main():
    if len(sys.argv) != 2:
        print('usage {0:s} fasta'.format(sys.argv[0]))
        sys.exit(1)
    fasta = sys.argv[1]
    seq1, seq2 = get_seq(fasta)
    ratio = count_transition(seq1, seq2) / count_transversion(seq1, seq2)
    print(ratio)
    sys.exit(0)


if __name__ == '__main__':
    main()
