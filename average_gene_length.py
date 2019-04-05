#! /usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--gff', required=True, help='gene annotation gff file')
args = parser.parse_args()


class GFFline:
    
    def __init__(self, gffline):
        l = gffline.split('\t')
        end = int(l[4])
        start = int(l[3])
        
        self.length = end - start
        if self.length < 1:
            print('calc wrong')
    
    
genelen = []
    
with open(args.gff, 'r') as infile:
    for line in infile:
        if line.startswith('NEMVE'):
            if line.split('\t')[2] == 'gene':
                g = GFFline(gffline=line)
                genelen.append(g.length)
            
            
average_gene_length = sum(genelen) / len(genelen)

print(average_gene_length)