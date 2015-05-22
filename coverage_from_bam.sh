#!/bin/bash

BED=$1
BAM=$2

# Get the length from the bed file
BED_LEN=$(awk '{if (/^[^#]/) {d[$1]+=($3-$2)}} END {for (c in d) {s+=d[c]}; print  s}' $BED)

# Get the total coverage from the bam file
BAM_COV=$(samtools mpileup -d 100000 -Q 5 -l $BED $BAM | cut -f 4 | datamash sum 1)
test -z "$BAM_COV" && BAM_COV=0

echo "BED length: $BED_LEN"
echo "BAM total: $BAM_COV"
AVG_COV=$(echo $BAM_COV/$BED_LEN | bc -l)
echo "Average coverage: $AVG_COV"
