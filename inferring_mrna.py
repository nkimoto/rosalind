#!usr/bin/python
import sys

def return_codon_num(rna):
    if rna in ('M', 'W'):
        return 1
    elif rna in ('F', 'Y', 'C', 'H', 'Q', 'N', 'K', 'D', 'E'):
        return 2
    elif rna == 'I':
        return 3
    elif rna in ('P', 'T', 'A', 'V', 'G'):
        return 4
    elif rna in ('L', 'R', 'S'):
        return 6
    else:
        print('It`s not mRNA!!')
        sys.exit()

def calculate_million_mod(num_list):
    res = 1
    for i in num_list:
        res = res * i
        while res >= 1000000:
            res = res%1000000
    return res
        

if __name__ == "__main__":
    rna_strings = sys.argv[1]
    num_list = []
    for i in rna_strings:
        num_list.append(return_codon_num(i))
    print(calculate_million_mod(num_list) * 3 % 1000000)
