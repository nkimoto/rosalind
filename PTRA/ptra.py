#!/usr/bin/env python

##############################################################################
#
# Name: needle.py
# Date: 2018-09-11
# Author: kimoton
# Description:
# Problem
# The 20 commonly occurring amino acids are abbreviated by using 20 letters from the English alphabet (all letters except for B, J, O, U, X, and Z). Protein strings are constructed from these 20 symbols. The RNA codon table shows the encoding from each RNA codon to the amino acid alphabet.
#
# The Translate tool from the SMS 2 package can be found here in the SMS 2 package
#
# A detailed list of genetic code variants (codon tables) along with indexes representing these codes (1 = standard genetic code, etc.) can be obtained here.
#
# For now, when translating DNA and RNA strings, we will start with the first letter of the string and ignore stop codons.
#
# Given: A DNA string s of length at most 10 kbp, and a protein string translated by s.
#
# Return: The index of the genetic code variant that was used for translation. (If multiple solutions exist, you may return any one.)
#
# Sample Dataset
# ATGGCCATGGCGCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA
# MAMAPRTEINSTRING
# Sample Output
# 1
#
##############################################################################

import sys
from Bio.Seq import translate

TABLE_NUM = (1, 2, 3, 4, 5, 6, 9, 10, 11, 12,
             13, 14, 16, 21, 23, 24, 25, 26,
             27, 28, 29, 30, 31)


def getTrueTable(seq, protein_seq):
    protein_seq = protein_seq.strip() + '*'
    print('Ans : {}'.format(protein_seq))
    for i in TABLE_NUM:
        protein_query = translate(seq.strip(), table=i)
        if protein_seq == protein_query:
            return i
        else:
            print(protein_query)
            continue


def main():
    input_seq = sys.argv[1]
    with open(input_seq, 'r') as rf:
        lines = rf.readlines()
        dna_seq = lines[0]
        protein_seq = lines[1]
    table_id = getTrueTable(dna_seq, protein_seq)
    print(table_id)


if __name__ == '__main__':
    main()
