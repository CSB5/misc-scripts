#!/usr/bin/env python
"""
Convert FASTA to FASTQ file with a static
 
Usage:
$ ./fasta_to_fastq NAME.fasta NAME.fastq


Taken from 
https://www.biostars.org/p/99886/
https://gist.github.com/mdshw5/c7cf7a232b27de0d4b31#file-fasta_to_fastq-py
"""
 
import sys, os
from Bio import SeqIO
 
# Get inputs
try:
    fa_path = sys.argv[1]
    fq_path = sys.argv[2]
    q = int(sys.argv[3])
except IndexError:
    sys.stderr.write("usage: fa fq q\n")
    sys.exit(1)
    
# make fastq
if fq_path=="-":
    fastq = sys.stdout
else:
    fastq = open(fq_path, "w")

if fa_path=="-":
    fasta = sys.stdin
else:
    fasta = open(fa_path, "r")

for record in SeqIO.parse(fasta, "fasta"):
    record.letter_annotations["phred_quality"] = [q] * len(record)
    SeqIO.write(record, fastq, "fastq")
