#!/bin/bash

BED=$1
BAM=$2
test -z "$BED" && exit 1
test -z "$BAM" && exit 1
test -e "$BED" || exit 1
test -e "$BAM" || exit 1

# Get the length from the bed file
BED_LEN=$(awk '{if (/^[^#]/) {d[$1]+=($3-$2)}} END {for (c in d) {s+=d[c]}; print  s}' $BED)

# Get the total coverage from the bam file
BAM_COV=$(grep -v '^#' $BED | while read c s e _; do let s=s+1; samtools mpileup -d 100000 -Q 5 -r "$c:$s-$e" $BAM | cut -f 4; done | datamash sum 1)
test -z "$BAM_COV" && BAM_COV=0

echo "BED length: $BED_LEN"
echo "BAM total: $BAM_COV"
AVG_COV=$(echo $BAM_COV/$BED_LEN | bc -l)
echo "Average coverage: $AVG_COV"
