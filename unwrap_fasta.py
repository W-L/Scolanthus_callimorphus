#!/usr/bin/env python3

# transforms a multi-line fasta to single-line

import sys

def unwrap_fasta(fa):
    f = open(fa, 'r')
    fu_name = fa + ".unwrapped"
    fu = open(fu_name, 'w+')
    seq = ''

    for line in f:
        # when line is a header, write last seq and new header
        if line[0] == ">":
            if len(seq) is not 0:
                seq = seq.replace('\n', '') + '\n'
                fu.write(seq)
            seq = ''
            fu.write(line)
        else:
            seq = seq + line

    # capture last seq when no more headers in file
    seq = seq.replace('\n', '') + '\n'
    fu.write(seq)

    f.close()
    fu.close()


for i in sys.argv[1:]:
    unwrap_fasta(i)
