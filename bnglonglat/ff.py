#!/usr/bin/env python3
#
# Pure Python implementation of the convertbng.util.convert_lonlat method.
#
with open('ostn15.rs', 'r') as f:

    n = n1 = n2 = 0
    pre = True

    for i, line in enumerate(f):
        n += 1
        if line.startswith('        ('):
            if pre:
                n1 += 1
            else:
                n2 += 1
        else:
            print(f'{i}: {line.strip()}')
            if 'entries' in line:
                pre = False
    print(n1, n2)
    print(n2 / n1)
