#! /usr/bin/env python3

import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('--reapr_out', required=True, help='reapr output called 01.stats.per_base.gz')
parser.add_argument('--feature', required=True, help='feature to convert to WIG')
args = parser.parse_args()


class REAPRline():
    
    'chr    pos     perfect_cov     read_cov        prop_cov        orphan_cov \
    bad_insert_cov  bad_orient_cov  read_cov_r      prop_cov_r      orphan_cov_r \
    bad_insert_cov_r        bad_orient_cov_r        frag_cov  frag_cov_err    FCD_mean  \
    clip_fl clip_rl clip_fr clip_rr FCD_err mean_frag_length'
    
    def __init__(self, reapr_line):
        l = reapr_line.rstrip('\n').split('\t')

        self.chr, self.pos, self.perfect_cov, self.read_cov, self.prop_cov = l[0:5]
        self.orphan_cov, self.bad_insert_cov, self.bad_orient_cov, self.read_cov_r, self.prop_cov_r = l[5:10]
        self.orphan_cov_r, self.bad_insert_cov_r, self.bad_orient_cov_r, self.frag_cov, self.frag_cov_err = l[10:15]
        self.FCD_mean, self.clip_fl, self.clip_rl, self.clip_fr, self.clip_rr = l[15:20]
        self.FCD_err, self.mean_frag_length = l[20:23]
        
        self.chr = self.chr.split('_')[0]
        
    
with gzip.open(args.reapr_out, 'rt') as infile:
    
    header = infile.readline()
    # write the track type and name
    # line = f"track type=wiggle_0 name={args.reapr_out}.{args.feature}\n"
    
        
    contig = "init"
    for line in infile:
        curr_line = REAPRline(reapr_line=line)
            
        # if the current chr changes -> new section of WIG file
        if curr_line.chr != contig:
            contig = curr_line.chr
            outline = f"variableStep chrom={contig}"
            print(outline)
                
        val = f"{curr_line.pos}\t{getattr(curr_line, args.feature)}"
        print(val)
            
            
            
                
            
            
            
        
        