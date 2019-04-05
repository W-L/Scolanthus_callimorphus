#! /usr/bin/env python3

import argparse
from glob import glob
from collections import defaultdict
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--aug_out', required=True, help='folder with augustus output')
args = parser.parse_args()


class Augfile:
    
    def __init__(self, name):
        self.name = name
        self.contig = name.split('/')[1].split('.')[0]
        self.genes = []
        
    def read_genes(self):
        genenum = 0
        raw_gene = ''
        with open(self.name, 'r') as infile:
            for line in infile:
                if line.startswith('# start gene '):
                    genenum = line.rstrip('\n').split('# start gene ')[1][1:]
                    
                    nline = infile.readline()

                    while nline.startswith('# end gene') is False:
                        raw_gene += nline
                        nline = infile.readline()
                        

                    gene = Gene(raw=raw_gene, contig=self.contig, num=genenum)
                    
                    self.genes.append(gene)
                    # reset the gene content
                    raw_gene = ''
                
                    

    
class Gene:
    
    def __init__(self, raw, contig, num):
        self.contig = contig
        self.num = num
        
        raw_split = raw.split('protein sequence = [')[1]
        raw_split2 = raw_split.split(']\n# Evidence')[0]
        prot_seq = raw_split2.replace('\n# ', '')

        self.seq = prot_seq
        
    def make_header(self):
        head = f'>{self.contig}_{self.num}'
        return(head)
    

    
aug_files = glob(args.aug_out + '/*.aug')


sequences = defaultdict(str)

for f in aug_files:
    augfile = Augfile(name=f)
    augfile.read_genes()
    
    for g in augfile.genes:
        g_head = g.make_header()
        if g_head not in sequences.keys():
            sequences[g_head] = g.seq
        else:
            print('sequence head is not unique' + g_head)
            sys.exit()

for header, seq in sequences.items():
    print(header)
    print(seq)
    


