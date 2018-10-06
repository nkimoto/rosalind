#!/usr/bin/env python

##############################################################################
#
# Name: ba9h.py
# Date: 2018-10-06
# Author: kimoton
# Description:
# Problem
# A version of FastQC can be downloaded here and run locally on any operating
# system with a suitable Java Runtime Environment (JRE) installed.

# An online version of FastQC is also available here in the "Andromeda" Galaxy
# instance.

# Given: A quality threshold, along with FASTQ entries for multiple reads.

# Return: The number of reads whose average quality is below the threshold.

# Sample Dataset
# 28
# @Rosalind_0041
# GGCCGGTCTATTTACGTTCTCACCCGACGTGACGTACGGTCC
# +
# 6.3536354;.151<211/0?::6/-2051)-*"40/.,+%)
# @Rosalind_0041
# TCGTATGCGTAGCACTTGGTACAGGAAGTGAACATCCAGGAT
# +
# AH@FGGGJ<GB<<9:GD=D@GG9=?A@DC=;:?>839/4856
# @Rosalind_0041
# ATTCGGTAATTGGCGTGAATCTGTTCTGACTGATAGAGACAA
# +
# @DJEJEA?JHJ@8?F?IA3=;8@C95=;=?;>D/:;74792.
# Sample Output
# 1
#
##############################################################################

import sys
from Bio import SeqIO


def read_avg_qual(read_record):
    qual_list = read_record.letter_annotations['phred_quality']
    return sum(qual_list) / len(qual_list)


def main():
    fastq = sys.argv[1]
    threshold = float(sys.argv[2])
    count = 0
    for read in SeqIO.parse(fastq, "fastq"):
        avg_qual = read_avg_qual(read)
        if avg_qual < threshold:
            count += 1
    return count


if __name__ == '__main__':
    print(main())
