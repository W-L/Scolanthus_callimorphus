#! /usr/bin/env python3

import argparse
from glob import glob
import ParseMyAugustus

parser = argparse.ArgumentParser()
parser.add_argument('--aug_out', required=True, help='folder with augustus output')
args = parser.parse_args()


# glob all augustus files in a directory
aug_files = glob(args.aug_out + '/*.aug')


for f in aug_files:
    # load the augustusfile
    augfile = ParseMyAugustus.Augfile(name=f)
    augfile.read_genes()

    for g in augfile.genes:
        # grab the features of each gene to concat the gtfs
        g.grabFeatures()
        g.writeFeatures()
        


