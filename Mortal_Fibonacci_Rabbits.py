#!/usr/bin/env python
# -*- coding: utf-8 -*-


def input_file(file_name):
    with open(file_name,'r') as rf:
        conts = [k.strip() for k in rf.read().split(' ')]
    return conts


def fib(n,m):
    """ 
    n : after n-months
    m : die after m-months
    """
    #initialize
    generation  = [0]*int(m)
    generation[0],generation[1] = 0,1

    for i in range(2,int(n)):
        temp = list(generation)
        generation[0] = sum(generation[1:])
        for j in range(1,len(generation)):
            generation[j] = temp[j-1]
        print(temp)
    return sum(generation)

if __name__ == "__main__":
    n,m = input_file("rosalind_fibd.txt")
    print(fib(n,m))

