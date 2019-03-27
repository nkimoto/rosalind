#!/usr/bin/env julia

#=
#
# Name: lexf.jl
# Date: 2019-03-27
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
=#

using ArgParse
using Base

"""
get all possible permutations by components_list
"""
function get_lines(input_file)
    rf = open(input_file, "r")
    lines = readlines(rf)
    return lines
end


"""
get all possible permutations by components_list
"""
function get_all_permutations(component_list, length)
    args = Iterators.repeat([component_list], length)
    all_perm = vec(collect(Iterators.product(args...)))
    sorted_all_perm = sort(all_perm)
    return sorted_all_perm
end


function main()
    s = ArgParseSettings()
    @add_arg_table s begin
        "input"
        help = "input file"
        required = true
    end
    args = parse_args(ARGS, s)
    input = args["input"]
    lines = get_lines(input)
    component_list = [String(i) for i in split(lines[1], " ")]
    length = parse(Int64, lines[2])
    permutations = get_all_permutations(component_list, length)
    println(join([join(i, "") for i in sort(permutations)], "\n"))
end


if occursin(PROGRAM_FILE, @__FILE__)
    main()
end
