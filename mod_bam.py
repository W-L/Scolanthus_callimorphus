#! /usr/bin/env python3

import argparse
import pysam


parser = argparse.ArgumentParser()
parser.add_argument('--bam', required=True, help='bam file')
args = parser.parse_args()


bamfile = pysam.AlignmentFile(args.bam, "rb")

# get header names and lengths
reflist = bamfile.header.references
lenlist = list(bamfile.header.lengths)

# change header names 
references = [r.split('_')[0] for r in reflist]

# construct new alignmentHeader object
new_header = pysam.AlignmentHeader.from_references(reference_names=references, reference_lengths=lenlist)

# write new bam
newbam = pysam.AlignmentFile(args.bam + '.new', "wb", header=new_header)

for read in bamfile:
    newbam.write(read)


newbam.close()
bamfile.close()