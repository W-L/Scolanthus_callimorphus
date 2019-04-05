#! /usr/bin/env python3

import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('--assembly', required=True, help='assembly file')
args = parser.parse_args()

outname = args.assembly.split('/')[-1]

with open(args.assembly, 'r') as infile:
    with open(f'{outname}.contig_sizes', 'w') as outfile:
        
        header = 'init'
        seqlen = 0
            
        for line in infile:
            if line.startswith('>'):
                
                # write previous contig
                outline = f'{header}\t{seqlen}\n'
                if header != 'init':
                    outfile.write(outline)
                
                # get new name, reset counter
                header = line.lstrip('>').rstrip('\n').split('_')[0]
                seqlen = 0
                continue
                    
            seqlen += len(line.rstrip('\n'))
                
        # write last contig
        outline = f'{header}\t{seqlen}\n'
        outfile.write(outline)

