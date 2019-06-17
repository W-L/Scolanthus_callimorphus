#! /usr/bin/env python3

import argparse
from collections import defaultdict
import subprocess
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--blastout', required=True, help='blast output from busco')
parser.add_argument('--contigs', required=True, help='folder with singulated contigs')
args = parser.parse_args()


def find_redundancy(blastout):
	# takes the blastfile output from a busco run and returns a dict of
	# possibly redundant contigs

	hits = defaultdict(list)

	with open(blastout, 'r') as blastfile:
		for line in blastfile:
			if line.startswith('#'):
				continue
			else:
				l = line.split('\t')
				if l[1] not in hits[l[0]]:
					hits[l[0]].append(l[1])
	return(hits)


def mummer_contigs(contigs, contig1, contig2):
	# takes two contigs and runs mummer
	name = contig1 + '_' + contig2
	
	c1 = f"{contigs}/{contig1}.fa"
	c2 = f"{contigs}/{contig2}.fa"
	
	task = ["../../builds/MUMmer3.23/nucmer", "--maxmatch", "-p", name, c1, c2]
	
	running = subprocess.Popen(' '.join(task), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
					 encoding='utf-8', shell=True)
	stdout, stderr = running.communicate()
	
	print(stdout)
	print(stderr)
	
	task = ["../../builds/MUMmer3.23/mummerplot", "-layout", "-large", "-filter", "--png", "-p", f"{name}.plot", f"{name}.delta"]
	
	running = subprocess.Popen(' '.join(task), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
					 encoding='utf-8', shell=True)
	stdout, stderr = running.communicate()
	
	print(stdout)
	print(stderr)

reduntigs = find_redundancy(blastout=args.blastout)
# check number of items in reduntigs elements

for gene, clist in reduntigs.items():
	print(len(clist))
	if len(clist) == 2:
		mummer_contigs(contigs=args.contigs, contig1=clist[0], contig2=clist[1])
	
	