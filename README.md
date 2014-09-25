misc-scripts
============

A collection of miscellaneous scripts used in CSB5

Script                     | Description
------                     |  -----------
bam_num_reads.sh           | Infer number of aligned and unaligned reads from BAM file[s]
bwamem_pe.sh               | Perform BWA-MEM mapping for paired end reads (incl. some massaging)
deinterleave_fastq.py      | Deinterleave an interleaved FastQ file
dwgsim_correctly_mapped.py | Extract only correctly mapped dwgsim reads from BAM file
fastq_num_reads.sh         | Infer number of reads from FastQ file[s]
ngs_pipeline.py            | Run markduplicates, indel realignment and base-quality recalibration
samflag.py                 | Explain SAM bit flags
seqgrep.py                 | Grep for sequence files searching name or sequence
fasta_to_fastq.py          | Convert fasta to fastq file with uniformly set quality

If not explicitely mentioned in the source file, the license is GPL2. See LICENSE
