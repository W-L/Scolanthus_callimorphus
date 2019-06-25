#! /usr/bin/env bash

genome=$1
index=$2
conserved_reads1=/scratch/weilguny/processing/16_readsFromConservedGenes/conservedReads_1.fastq
conserved_reads2=/scratch/weilguny/processing/16_readsFromConservedGenes/conservedReads_2.fastq
conserved_genes=/scratch/weilguny/processing/16_readsFromConservedGenes/scol_trinity_conserved_genes.fasta

# 1 map reads of conserved genes to the assembly
# make a genome index first
mkdir -p $index
/proj/rpz/slurm_scripts/tmpfile_env --in $genome "STAR --runMode genomeGenerate --runThreadN 8 --genomeDir $index --genomeFastaFiles %i1" >STAR_genome_gen.job

# then map reduced read set to assembly
/proj/rpz/slurm_scripts/tmpfile_env --in $index --in $conserved_reads1 --in $conserved_reads2 "STAR --runThreadN 8 --outSAMtype BAM Unsorted SortedByCoordinate --genomeDir %i1 --readFilesIn %i2 %i3 --outFileNamePrefix ${index}_star_ --outFilterType BySJout --outFilterMultimapNmax 100 --alignSJoverhangMin 8 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 500000 --alignMatesGapMax 500000 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3" >STAR_map.job

# 2 map conserved gene seqs from trinity assebmly to genome assembly
# build genome for gmap first
/proj/rpz/slurm_scripts/tmpfile_env --in $genome "gmap_build -D ${index}_GMAP -d ${index}_GMAP %i1" >gmap_build.job

# then map
/proj/rpz/slurm_scripts/tmpfile_env --in ${index}_GMAP --in $conserved_genes "gmap -D %i1 -d ${index}_GMAP --split-large-introns --min-intronlength=20 --min-identity=.97 --min-trimmed-coverage=.8 --npaths=1 --no-chimeras --format psl -t 8 %i2 > ${index}_conserved_genes.psl" >gmap_genes.job

