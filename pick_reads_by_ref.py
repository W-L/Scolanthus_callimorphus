#! /usr/bin/env python3

import argparse
import pysam


parser = argparse.ArgumentParser()
parser.add_argument('--bam', required=True, help='bam file')
parser.add_argument('--reflist', required=True, help='list of references from which to extract reads')
args = parser.parse_args()


# read the list of references
refs = []
with open(args.reflist, 'r') as reflist:
    for line in reflist:
        refs.append(line.rstrip('\n'))
print(f'found {len(refs)} references')

# open the original bam file
bamfile = pysam.AlignmentFile(args.bam, "rb")

# open new bam, with same header (template=)
newbam = pysam.AlignmentFile(args.bam + '.extracted', "wb", template=bamfile)

# for every read in original bam - check if aligned to one of the refs
written=0
skipped=0
for read in bamfile:
    if read.reference_name in refs:
        newbam.write(read)
        written += 1
    else:
        skipped += 1

print(f'wrote {written} reads')
print(f'skipped {skipped} reads')

newbam.close()
bamfile.close()