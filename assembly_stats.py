#! /usr/bin/env python3

import argparse
import numpy as np
from collections import defaultdict


parser = argparse.ArgumentParser()
parser.add_argument('--assemblies', required=True, help='comma separated assemblies')
parser.add_argument('--genomesize', help='estimate of the genome size for normalization', default=400000000)
parser.add_argument('--output', required=True, help='name for output file')
parser.add_argument('--avgenelen', help='average gene length', default=4800)
args = parser.parse_args()


class Assembly:
    
    def __init__(self, name):
        self.name = a
        
    def count_contig_lengths(self):
        self.contigs = []
        with open(self.name, 'r') as infile:
            seqlen = 0
            
            for line in infile:
                if line.startswith('>'):
                    self.contigs.append(seqlen)
                    seqlen = 0
                    continue
                
                seqlen += len(line.rstrip('\n'))
        
        self.contigs.sort()
        
        # make another copy with reverse sorting for cumsum
        self.rev_contigs = list(self.contigs)
        self.rev_contigs.sort(reverse=True)
        
    def contig_cumsum(self):
        cumsum = list(np.cumsum(self.rev_contigs))
        return(cumsum)
        
    
    def calcNG(self, genomesize):
        self.NGs = []
        
        for i in range(1, 101):
            contigs = self.contigs.copy()
            threshold = genomesize * (i / 100)
            seq_sum = 0
            
            while seq_sum < threshold:
                if len(contigs) != 0:
                    curr_seq = contigs.pop()
                    seq_sum += curr_seq
                    
                else:
                    curr_seq = 0
                    break
            
            self.NGs.append(curr_seq)
            
    def perc_genome_bigger_than_genes(self, genomesize, average_gene_length=1):
        contigs = self.contigs.copy()
        
        big_contigs = [c for c in contigs if c > average_gene_length]
        
        percentage = sum(big_contigs) / genomesize
        self.percentage = percentage * 100
        
        
def write_NGs(NGdict, out):
    with open(out, 'w') as outfile:
        keys = list(NGdict.keys())
        header = 'NGx ' + (' ').join(keys) + '\n'
        outfile.write(header)
        
        lines = zip(*list(NGdict.values()))
        c = 1
        for l in lines:
            joined_line = (' ').join([str(x) for x in l])
            row = str(c) + ' ' + joined_line + '\n'
            c += 1
            outfile.write(row)
            
            
def write_dict(dicti, out):
    with open(out, 'w') as outfile:
        for key, val in dicti.items():
            line = str(key) + ' ' + str(val) + '\n'
            outfile.write(line)
            
def write_vals(dicti, out):
    with open(out, 'w') as outfile:
        for key, val in dicti.items():
            line = str(key) + ' ' + str(' '.join([str(i) for i in val])) + '\n'
            outfile.write(line)
            

assem = args.assemblies.split(',')

NG_vals = defaultdict(list)
percentages = defaultdict(int)
cumsums = defaultdict(int)

for a in assem:
    # read the assembly
    ass = (Assembly(name=a))
    # count the length of all contigs
    ass.count_contig_lengths()
    # transfer cumsum
    cumsums[ass.name] = ass.contig_cumsum()
    # calculate NG 1-100
    ass.calcNG(genomesize=int(args.genomesize))
    # transfer to a collection list
    NG_vals[ass.name] = ass.NGs
    
    # calc percentage of genome in contigs > average gene
    ass.perc_genome_bigger_than_genes(average_gene_length=int(args.avgenelen),
                                      genomesize=int(args.genomesize))             
    # collect in dict                                  
    percentages[ass.name] = ass.percentage
    
    
write_NGs(NGdict=NG_vals, out=args.output + '.NG')
write_dict(dicti=percentages, out=args.output + '.perc')
write_vals(dicti=cumsums, out=args.output + '.cumsum')

