#!/usr/bin/env python

##############################################################################
#
# Name: lexf.py
# Date: 2019-03-23
# Author: kimoton
# Description:
# Problem
# Assume that an alphabet A has a predetermined order; that is, we write the
# alphabet as a permutation A=(a1,a2,…,ak), where a1<a2<⋯<ak. For instance, the
# English alphabet is organized as (A,B,…,Z).
#
# Given two strings s and t having the same length n, we say that s precedes
# t in the lexicographic order (and write s<Lext) if the first symbol s[j]
# that doesn't match t[j] satisfies sj<tj in A.
#
# Given: A collection of at most 10 symbols defining an ordered alphabet, and a
# positive integer n (n≤10).
#
# Return: All strings of length n that can be formed from the alphabet, ordered
# lexicographically (use the standard order of symbols in the English
# alphabet).
#
# Sample Dataset
# A C G T
# 2
# Sample Output
# AA
# AC
# AG
# AT
# CA
# CC
# CG
# CT
# GA
# GC
# GG
# GT
# TA
# TC
# TG
# TT
##############################################################################

import sys
import itertools
import operator


def get_lines(_file):
    with open(_file, "r") as rf:
        lines = rf.readlines()
    return lines


def get_all_permutations(components_list, length):
    """
    get all possible permutations by components_list
    """
    all_perm = list(itertools.product(components_list, repeat=length))
    sorted_all_perm = sorted(all_perm)
    return sorted_all_perm


def main():
    if len(sys.argv) != 2:
        print('usage {0:s} input_file'.format(sys.argv[0]))
        sys.exit(1)
    input_file = sys.argv[1]
    line1, line2 = get_lines(input_file)
    components_list = line1.split()
    length = int(line2)
    permutations = get_all_permutations(components_list, length)
    for i in permutations:
        print("".join(i))
    sys.exit(0)


if __name__ == '__main__':
    main()
