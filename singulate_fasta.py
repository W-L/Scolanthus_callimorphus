#!/usr/bin/env python3

# save fasta seqs individually

import sys

def singulate_fasta(fa):
    with open(fa, 'r') as fastafile:
        for line in fastafile:
            if line.startswith('>'):
                file = line.lstrip('>').rstrip('\n')
                outfile = open(file + '.fa', 'w')
                outfile.write(line)
                continue
            outfile.write(line)
            outfile.close()
    

singulate_fasta(sys.argv[1])

