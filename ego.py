#! usr/bin/env python

import itertools

perms = [str(i)[1:-1] for i in itertools.permutations(range(1,7))]
print(len(perms))
for i in perms:
    print(''.join(i.split(',')))



