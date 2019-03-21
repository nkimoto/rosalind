#!/usr/bin/env python

##############################################################################
#
# Name: needle.py
# Date: 2019-03-21
# Author: kimoton
# Description:
# Problem
# After identifying the exons and introns of an RNA string, we only need to
# delete the introns and concatenate the exons to form a new string ready for
# translation.
#
# Given: A DNA string s (of length at most 1 kbp) and a collection of
# substrings of s acting as introns. All strings are given in FASTA format.
#
# Return: A protein string resulting from transcribing and translating the
# exons of s. (Note: Only one solution will exist for the dataset provided.)
#
# Sample Dataset
# >Rosalind_10
# ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG
# >Rosalind_12
# ATCGGTCGAA
# >Rosalind_15
# ATCGGTCGAGCGTGT
# Sample Output
# MVYIADKQHVASREAYGHMFKVCA
#
##############################################################################

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_rna

def get_seq(fasta):
    records = []
    with open(fasta, "r") as rf:
        for record in SeqIO.parse(rf, "fasta"):
            records.append(str(record.seq))
    return records


def remove_intron(seq, intron_set):
    for intron in intron_set:
        seq = seq.replace(intron, "")
        print(seq)
    return seq


def main():
    if len(sys.argv) != 2:
        print('usage {0:s} fasta'.format(sys.argv[0]))
        sys.exit(1)
    fasta = sys.argv[1]
    seqs = get_seq(fasta)
    seq = remove_intron(seqs[0], seqs[1:])
    mrna = Seq(seq, generic_rna).translate()
    print(mrna)
    sys.exit(0)


if __name__ == '__main__':
    main()
