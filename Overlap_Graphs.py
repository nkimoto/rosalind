#!usr/bin/env python
# -*- coding: utf-8 -*-

import itertools
from Bio import SeqIO
import sys
from collections import OrderedDict


def read_file(rf):
    with open(rf,'r') as rf2:
        return rf2.read()

def read_seq(fasta):
    inp = SeqIO.parse(fasta,'fasta')
    pre_suf_dict = OrderedDict()
    edges = []
    for i in inp:
        pre_suf_dict[i.id] = [str(i.seq[:3]), str(i.seq[-3:])]
    for n,j in enumerate(pre_suf_dict):
        for k in list(pre_suf_dict.keys())[n+1:]:
            if pre_suf_dict[j][0] == pre_suf_dict[k][1] or pre_suf_dict[j][1] == pre_suf_dict[k][0]:
                edges.append((j,k))
    for l in edges:
        print(" ".join(l))

def read_seq2(fasta):
    inp = SeqIO.parse(fasta,"fasta")
    id_seq_dict = {}
    for i in inp:
        id_seq_dict[str(i.id)] = str(i.seq)
    return id_seq_dict

def is_k_overlap(s1,s2,k):
    return s1[-k:] == s2[:k]

def k_edges(data,k):
    edges = []
    for u,v in itertools.combinations(data,2):
        u_dna,v_dna = data[u], data[v]

        if is_k_overlap(u_dna, v_dna,k):
            edges.append((u,v))
        if is_k_overlap(v_dna, u_dna,k):
            edges.append((v,u))
    return edges


if __name__ == '__main__':
    args = sys.argv
    fasta = read_seq(args[1])
    print()
    fasta2 = read_seq2(args[1])
    for j in k_edges(fasta2,3):
        print(" ".join(j))
