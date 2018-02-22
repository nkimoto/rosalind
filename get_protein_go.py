from Bio import ExPASy
from Bio import SwissProt
input_ids = input().split('\n')
for i in input_ids:
    handle = ExPASy.get_sprot_raw(i)
    record = SwissProt.read(handle)
    print(dir(record))
    print(record.cross_references)
    print('\n'.join([r[2].split(':')[1] for r in record.cross_references if r[0] == 'GO' and r[2].split(':')[0] == 'P']))
