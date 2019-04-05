#! /usr/bin/env bash

module load augustus
export AUGUSTUS_CONFIG_PATH=./config

# set up the output
mkdir -p augustus_out

# declare folder with the single-seq fasta, no slash after name
splitFasta=$1

# get the contig files
contigs=$(ls -1 $splitFasta/*)

# for each contig
for c in $contigs;
do
	c_name=$(basename $c)
	augustus_cmd="augustus \
		--/Constant/min_coding_len=101 \
		--alternatives-from-evidence=true \
		--sample=0 \
		--species=sc1_canu_full \
		--extrinsicCfgFile=extrinsic.cfg \
		--hintsfile=sj.rm.hints \
		--UTR=off \
		$c \
		>augustus_out/$c_name.aug
"

	echo $augustus_cmd

done