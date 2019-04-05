#! /usr/bin/env bash

ill1=$1
ill2=$2

name=$(basename -s .ill1.bam $ill1)

mod="module load samtools"
eval $mod


head1="samtools view -H $ill1 >tmp.ill1.head"
head2="samtools view -H $ill2 >tmp.ill2.head"
# echo $head1
# echo $head2

eval $head1
eval $head2

body1="samtools view $ill1 >tmp.ill1.body"
body2="samtools view $ill2 >tmp.ill2.body"
# echo $body1
# echo $body2

eval $body1
eval $body2

diffhead="diff tmp.ill1.head tmp.ill2.head >tmp.diff"
# echo $diffhead
eval $diffhead

if [ -s tmp.diff ]
then
	echo "headers are not equal"
	exit 0
else
	# merge the two bam files with one header
	merge="cat tmp.ill1.head tmp.ill1.body tmp.ill2.body >$name.sam"
	# echo $merge
	eval $merge

	# sort the bamfile
	sorting="samtools view -h -b $name.sam | samtools sort -@ 4 -o $name.sort.bam -"
	# echo $sorting
	eval $sorting

	index="samtools index $name.sort.bam"
	# echo $index
	eval $index

fi