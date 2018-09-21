#!/usr/bin/env python

##############################################################################
#
# Name: needle.py
# Date: 2018-09-11
# Author: kimoton
# Description:
# Problem
# Problem
# An online interface to EMBOSS's Needle tool for aligning DNA and RNA strings can be found here.
#
# Use:
#
# The DNAfull scoring matrix; note that DNAfull uses IUPAC notation for ambiguous nucleotides.
# Gap opening penalty of 10.
# Gap extension penalty of 1.
# For our purposes, the "pair" output format will work fine; this format shows the two strings aligned at the bottom of the output file beneath some statistics about the alignment.

# Given: Two GenBank IDs.

# Return: The maximum global alignment score between the DNA strings associated with these IDs.

# Sample Dataset
# JX205496.1 JX469991.1
# Sample Output
# 257
#
##############################################################################

import sys
from operator import add


mono_mass_table = {
    'A': 71.03711,
    'C': 103.00919,
    'D': 115.02694,
    'E': 129.04259,
    'F': 147.06841,
    'G': 57.02146,
    'H': 137.05891,
    'I': 113.08406,
    'K': 128.09496,
    'L': 113.08406,
    'M': 131.04049,
    'N': 114.04293,
    'P': 97.05276,
    'Q': 128.05858,
    'R': 156.10111,
    'S': 87.03203,
    'T': 101.04768,
    'V': 99.06841,
    'W': 186.07931,
    'Y': 163.06333,
}


def get_weight(protein_seq_file):
    with open(protein_seq_file, "r") as rf:
        seq = rf.read().strip()
    weights = [mono_mass_table[i] for i in list(seq)]
    total_weight = sum(weights)
    return total_weight


def main():
    if len(sys.argv) != 2:
        print('usage {0:s} protein_seq_file'.format(sys.argv[0]))
        sys.exit(1)
    protein_seq_file = sys.argv[1]
    weight = get_weight(protein_seq_file)
    print(weight)
    sys.exit(0)


if __name__ == '__main__':
    main()
