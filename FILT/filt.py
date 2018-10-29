#!/usr/bin/env python

##############################################################################
#
# Name: filt.py
# Date: 2018-10-29
# Author: kimoton
# Description:
# Problem
# Problem
# Poor-quality reads can be filtered out using the FASTQ Quality Filter tool from
# the FASTX toolkit. A command-line version of FASTX can be downloaded for Linux
# or MacOS from its website. An online interface for the FASTQ Quality Filter is
# also available here within the Galaxy web platform.
#
# Given: A quality threshold value q, percentage of bases p, and set of FASTQ
# entries.
#
# Return: Number of reads in filtered FASTQ entries
#
# Sample Dataset
# 20 90
# @Rosalind_0049_1
# GCAGAGACCAGTAGATGTGTTTGCGGACGGTCGGGCTCCATGTGACACAG
# +
# FD@@;C<AI?4BA:=>C<G=:AE=><A??>764A8B797@A:58:527+,
# @Rosalind_0049_2
# AATGGGGGGGGGAGACAAAATACGGCTAAGGCAGGGGTCCTTGATGTCAT
# +
# 1<<65:793967<4:92568-34:.>1;2752)24')*15;1,.3*3+*!
# @Rosalind_0049_3
# ACCCCATACGGCGAGCGTCAGCATCTGATATCCTCTTTCAATCCTAGCTA
# +
# B:EI>JDB5=>DA?E6B@@CA?C;=;@@C:6D:3=@49;@87;::;;?8+
# Sample Output
# 2
#
##############################################################################

import os.path
import sys
import subprocess


def fastx_toolkit(q, p, i, o):
    """
    q : threshold value
    p : perecentage of bases
    i : input
    o : output
    """
    cmd = f'fastq_quality_filter -i {i} -o {o} -Q33 \
            -p {p} -q {q}'
    result = subprocess.call(cmd, shell=True)
    return result


def count_read_num(input_f):
    with open(input_f, 'r') as rf:
        read_num = len(rf.readlines())/4
    return int(read_num)


def main():
    input_f = sys.argv[1]
    quality = sys.argv[2]
    percentage = sys.argv[3]
    output_f = os.path.splitext(input_f)[0] + '_filtered.fq'
    fastx_toolkit(q=quality, p=percentage,
                  i=input_f, o=output_f)
    read_num = count_read_num(output_f)
    print(read_num)
    sys.exit(0)


if __name__ == '__main__':
    main()
