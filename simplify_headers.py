#! /usr/bin/env python3

import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--fa', required=True, help='fasta file')
parser.add_argument('--delimiter', required=True, help='delimiter in the header')
args = parser.parse_args()


with open(args.fa, 'r') as infile:
    for line in infile:
        if line.startswith('>'):
            header = line.split(args.delimiter)[0]
            print(header)
            continue
            
        line = line.rstrip('\n')
        print(line)
        
        