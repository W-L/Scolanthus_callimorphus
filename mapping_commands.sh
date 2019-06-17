#! /usr/bin/env bash

load="module load star"
eval $load

genomes=$(ls *.fasta)

for g in $genomes; do

	MN=$(basename -s .fasta $g)
	
	mkdir -p $MN

	STAR --runMode genomeGenerate --runThreadN 8 --genomeDir $MN --genomeFastaFiles $g
	
done


for g in $genomes; do

	MN=$(basename -s .fasta $g)

	STAR --genomeDir $MN --readFilesIn $MN.fastq --outFileNamePrefix $MN.out --runThreadN 8

done
