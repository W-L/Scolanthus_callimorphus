#! /usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--fasta', required=True, help='fasta file')
parser.add_argument('--reflist', required=True, help='list of references to extract')
args = parser.parse_args()


refs = set()

# extract the reference names
with open(args.reflist) as reflist:
    for line in reflist:
        refs.add(line.rstrip('\n'))
        

# open the fasta file
with open(args.fasta) as fasta:
    for line in fasta:
        # for every header
        if line.startswith('>'):
            # check if header is in the refs
            header = line.split(' ')[0].lstrip('>')
            if header in refs:
                print(f'>{header}')
                l = fasta.readline()
                while l.startswith('>') is False:
                    
                    print(f'{l}', end='', flush=True)
                    l = fasta.readline()
            
                
            
        
        
# ../15_reducedRNAset/complete_buscos_from_scol_trans_sorted_dedup
# ../15_reducedRNAset/scol_trinity_reduced.fasta