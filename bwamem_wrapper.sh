#!/bin/bash

# This script will align your reads with BWA-MEM and massage the
# output nicely, by 1) adding a read-group, clean the output, fix
# the # mate information (if any) and sort the BAM file

# Authors:
# - Andreas Wilm <wilma@gis.a-star.edu.sg>
# - LI Chenhao <lich@gis.a-star.edu.sg>
#
# License: MIT

set -o pipefail

DEFAULT_THREADS=8

PICARDDIR_DEFAULT=/mnt/software/src/MAPPERS/picard-tools-1.113/
test -z "$PICARD_DIR" && PICARD_DIR=$PICARDDIR_DEFAULT
PICARD_ADDRG=${PICARD_DIR}/AddOrReplaceReadGroups.jar
PICARD_CLEAN=${PICARD_DIR}/CleanSam.jar
PICARD_FIXMATE=${PICARD_DIR}/FixMateInformation.jar

usage() {
cat <<EOF
$(basename $0): wrapper for running BWA-MEM on paired-end reads
    -f | --ref     : reference
    -1 | --fq1     : first fastq file
    -2 | --fq2     : second fastq file (optional)
    -o | --outpref : output prefix
    -t | --threads : number of threads to use (default is $DEFAULT_THREADS)
EOF
}


# check for programs
for b in bwa java; do
    if ! which $b >/dev/null 2>&1; then 
        echo "FATAL: $b seems to be missing on your system" 1>&2
        exit 1
    fi
done
# make sure bwa supportst mem
if bwa mem 2>&1 | grep -q unrecognized; then
    echo "FATAL: bwa doesn't seem to support the 'mem' command" 1>&2
    exit 1
fi
# check for picard stuff
for f in $PICARD_ADDRG $PICARD_CLEAN $PICARD_FIXMATE; do
    if [ ! -e "$f" ]; then
        echo "FATAL: missing file $f" 1>&2
        exit 1
    fi
done


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
JAVA_EXTRA_ARGS="-XX:ParallelGCThreads=$threads -Xmx4g"
# smetimes needed see https://sourceforge.net/p/samtools/mailman/message/31865651/
#JAVA_EXTRA_ARGS="$JAVA_EXTRA_ARGS -Dsamjdk.try_use_intel_deflater=false"
for f in $fq1 $ref; do
    test -z "$f" && exit 1
    test -e "$f" || exit 1
done

test -z "$outpref" && exit 1
test -e ${OUTPREF}.bam && exit 1



if [ ! -e ${ref}.bwt ]; then
    # childish attempt to avoid race condition
    sleep=$RANDOM
    let "sleep %= 10"
    sleep $sleep
    bwa index $ref || exit 1
fi


# use md5sum of the file name as the read group
rgid=$(echo $PWD/$fq1 | md5sum | cut -d " " -f1)
rgpu=$rgid.PU

PICARD_TMP=$(dirname $outpref)
# using /dev/shm/ only makes sense if RAM is big enough
#if test -d /dev/shm/; then
#	PICARD_TMP=/dev/shm/$PICARD_TMP
#fi
#cat <<EOF
bwa mem -M -t $threads $ref $fq1 $fq2 \
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
samtools index ${outpref}.bam
#EOF
