#!/bin/bash

set -o pipefail

# Author: Andreas Wilm <wilma@gis.a-star.edu.sg>
# License: WTFPL http://www.wtfpl.net/

# Determine number of reads contained in given (indexed) BAM files

echo -e "#file\t#aln\t#unaln\t#total"
for bam in "$@"; do
    echo -ne "$bam\t"
    samtools idxstats $bam | awk '{a+=$3; u+=$4} END {printf "%d\t%d\t%d\n", a, u, a+u}';
done
