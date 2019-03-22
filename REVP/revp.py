#!/usr/bin/env python

##############################################################################
#
# Name: tran.py
# Date: 2019-03-22
# Author: kimoton
# Description:
# Problem
#
# Figure 2. Palindromic recognition site
# A DNA string is a reverse palindrome if it is equal to its reverse
# complement. For instance, GCATGC is a reverse palindrome because its reverse
# complement is GCATGC. See Figure 2.
#
# Given: A DNA string of length at most 1 kbp in FASTA format.
#
# Return: The position and length of every reverse palindrome in the string
# having length between 4 and 12. You may return these pairs in any order.
#
# Sample Dataset
# >Rosalind_24
# TCAATGCATGCGGGTCTATATGCAT
# Sample Output
# 4 6
# 5 4
# 6 6
# 7 4
# 17 4
# 18 4
# 20 6
# 21 4
#
##############################################################################

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def get_rev(seq):
    seq = Seq(seq, generic_dna)
    return str(seq.reverse_complement())


def get_seq(fasta):
    records = []
    with open(fasta, "r") as rf:
        for record in SeqIO.parse(rf, "fasta"):
            records.append(str(record.seq))
    return records


def isParindrome(seq):
    if seq == get_rev(seq):
        return True


def find_palindrome(seq, _min, _max):
    len_list = []
    for l in range(_min, _max + 1):
        for i in range(0, len(seq) - l + 1):
            if isParindrome(seq[i:i+l]):
                len_list.append((i+1, l))
    len_list.sort()
    return len_list


def main():
    if len(sys.argv) != 2:
        print('usage {0:s} fasta'.format(sys.argv[0]))
        sys.exit(1)
    fasta = sys.argv[1]
    seq1 = get_seq(fasta)[0]
    len_list = find_palindrome(seq1, 4, 12)
    for i, j in len_list:
        print(i, j)
    sys.exit(0)


if __name__ == '__main__':
    main()
