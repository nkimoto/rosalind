#!/usr/bin/env python

##############################################################################
#
# Name: needle.py
# Date: 2018-09-15
# Author: kimoton
# Description:
# Problem
# Sometimes it's necessary to convert data from FASTQ format to FASTA format. For
# example, you may want to perform a BLAST search using reads in FASTQ format
# obtained from your brand new Illumina Genome Analyzer.

# Links:
#     A FASTQ to FASTA converter can be accessed from the Sequence conversion
#     website

#     A free GUI converter developed by BlastStation is available here for
#     download or as an add-on to Google Chrome.

#     There is a FASTQ to FASTA converter in the Galaxy web platform. Note that
#     you should register in the Galaxy and upload your file prior to using this
#     tool.

#     Given: FASTQ file

#     Return: Corresponding FASTA records

#     Sample Dataset
#     @SEQ_ID
#     GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
#     +
#     !*((((***+))%%%++)(%%%%).1***-+*****))**55CCF>>>>>>CCCCCCC65
#     Sample Output
#     >SEQ_ID
#     GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
#
##############################################################################
import os
import sys
from Bio import SeqIO

def converter(fq):
    output_fa = os.path.splitext(fq)[0] + '.fa'
    SeqIO.convert(fq, 'fastq', output_fa, 'fasta')


def main():
    input_fq = sys.argv[1]
    converter(input_fq)


if __name__ == '__main__':
    main()
