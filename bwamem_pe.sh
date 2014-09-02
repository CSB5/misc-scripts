#!/bin/bash

# This script will align your paired-end reads with BWA-mem. It will
# also add a read-group (md5 of your filename), clean the SAM and fix
# mate information during which the BAM file will get sorted as well.

# Authors:
# - Andreas Wilm <wilma@gis.a-star.edu.sg>
# - LI Chenhao <lich@gis.a-star.edu.sg>
#
# License: GPL2

set -o pipefail

DEFAULT_THREADS=10

# FIXME make BWA user changeable
BWA=/mnt/software/stow/bwa-0.7.10/bin/bwa

JAVA_EXTRA_ARGS="-XX:ParallelGCThreads=$DEFAULT_THREADS -Xmx4g"

# FIXME make PICARD_DIR user changeable
PICARD_DIR=/mnt/software/src/MAPPERS/picard-tools-1.113/
PICARD_ADDRG=${PICARD_DIR}/AddOrReplaceReadGroups.jar
PICARD_CLEAN=${PICARD_DIR}/CleanSam.jar
PICARD_FIXMATE=${PICARD_DIR}/FixMateInformation.jar

usage() {
cat <<EOF
$(basename $0): bwa mem pe wrapper
    -f | --ref     : reference
    -1 | --fq1     : first fastq file
    -2 | --fq2     : second fastq file
    -o | --outpref : output prefix
EOF
}


while [ "$1" != "" ]; do
    case $1 in
        -f | --ref )
            shift
            ref=$1
            ;;
        -t | --threads )
            shift
            threads=$1
            ;;
        -1 | --fq1 )
            shift
            fq1=$1
            ;;
        -2 | --fq2 )
            shift
            fq2=$1
            ;;
        -o | --outpref )
            shift
            outpref=$1
            ;;
        * ) 
            echo "FATAL: unknown argument \"$1\""
            usage
            exit 1
    esac
    shift
done

test -z "$threads" && threads=$DEFAULT_THREADS
for f in $fq1 $fq2 $ref; do
    test -z "$f" && exit 1
    test -e "$f" || exit 1
done
test -z "$outpref" && exit 1
test -e ${OUTPREF}.bam && exit 1

for f in $PICARD_ADDRG $PICARD_CLEAN $PICARD_FIXMATE $BWA; do
    if [ ! -e "$f" ]; then
        echo "FATAL: missing file $f" 1>&2
        exit 1
    fi
done


if [ ! -e ${ref}.bwt ]; then
    # childish attempt to avoid race condition
    sleep=$RANDOM
    let "sleep %= 10"
    sleep $sleep
    $BWA index $ref || exit 1
fi


# use md5sum of the file name as the read group
rgid=$(echo $PWD/$fq1 | md5sum | cut -d" " -f1)
rgpu=$rgid.PU

PICARD_TMP=$(dirname $outpref)
#cat <<EOF
$BWA mem -M -t $threads $ref $fq1 $fq2 \
    -R "@RG\tID:${rgid}\tPL:illumina\tPU:${rgpu}\tLB:lb-dummy\tSM:sm-dummy"| \
    java $JAVA_EXTRA_ARGS -jar $PICARD_CLEAN \
    VALIDATION_STRINGENCY=LENIENT \
    COMPRESSION_LEVEL=0 \
    INPUT=/dev/stdin OUTPUT=/dev/stdout | \
    java $JAVA_EXTRA_ARGS -jar $PICARD_FIXMATE \
    VALIDATION_STRINGENCY=LENIENT \
    SORT_ORDER=coordinate \
    TMP_DIR=$PICARD_TMP \
    INPUT=/dev/stdin OUTPUT=${outpref}.bam || exit 1
#EOF

samtools index ${outpref}.bam
