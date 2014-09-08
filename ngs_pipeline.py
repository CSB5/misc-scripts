#!/usr/bin/env python
"""Ruffus based, loose implementation of the data pre-processing part
of the GATK best practices pipeline.

See
http://www.broadinstitute.org/gatk/guide/best-practices
and
http://www.ruffus.org.uk/

The pipeline is based on file-extensions. The current flow of things
is as follows:

0   .bam (indexed input bam file)
1.1 .mdup.bam
1.2 .mdup.bam.bai
2.1 .realn.intervals (using 1.1)
2.2 .realn.bam (using 1.1 and 2.1)
3.1 .realn.recal.table (using 2.2)
3.2 .realn.recal.bam (using 2.2 and 3.1)
"""


__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2014 Genome Institute of Singapore"
__license__ = "GPL2"


# Ruffus documentation:
# http://www.ruffus.org.uk/contents.html
#
# Some snippets cargo culted from:
# https://github.com/fsroque/NGS-pipeline/blob/master/pipeline_multisample.py
# and
# https://github.com/seandavi/ngs/blob/master/scripts/exomecapture.py
#
# NOTE: this was the first time I used ruffus. There is a good chance
# therefore that I misused it.


# --- standard library imports
#
import sys
import os
import argparse
import logging
#import json
import subprocess

#--- third-party imports
#
from ruffus import *

#--- project specific imports
#
# /


# global logger
#
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')

# FIXME: functions to load and save config
CFG = {
    'java': 'java',# path to java binary
    'java_opts': ['-Xmx4g'],# has to be list
    'gatk': '/mnt/software/src/GATK/GenomeAnalysisTK-2.7-4-g6f46d11/GenomeAnalysisTK.jar',# path to GATK jar
    'picard_dir': '/mnt/software/src/MAPPERS/picard-tools-1.113/',
    'numthreads': 8,# number of threads to use
    'samtools': 'samtools',# path to samtools
    'trim_galore': 'trim_galore',# path to trim_galore
    'fastqc': 'fastqc',# path to fastqc
    'bam_in': [],# list of BAM input files
    'reffa': None,# reference fasta file
    'dbsnp': None,# dbsnp vcf files
    'simul': False}# simulate only

DEFAULT_TASK = 'baserecalibrator'

class JobFailedException(Exception):
    pass


def run_cmd(cmd, log_stdout=sys.stdout, log_stderr=sys.stderr):
    """Wrapper to Popen. Will raise an error if command execution
    failed with non-zero exit code
    """

    p = subprocess.Popen(cmd,# shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    (stdout, stderr) = p.communicate()

    for line in stdout.splitlines():
        log_stdout.write("stdout: %s\n" % line)
    for line in stderr.splitlines():
        log_stderr.write("stderr: %s\n" % line)
    # FIXME will only print stdout and stderr once done and also keep
    # all the output in memory right? Better to use shell and redirect
    # instead?

    if p.returncode != 0:
        raise JobFailedException(
            "Following command failed with exit status"
            " %s: '%s'. stderr was '%s'" % (
                p.returncode, ' '.join(cmd), stderr))

    return (stdout.splitlines(), stderr.splitlines())

# ----------------------------------------------------------------------
# SO unused bits
# ----------------------------------------------------------------------

def run_trimming(fastq_files, phred33=True, outdir=None):
    """Run FastQ trimmer

    fastq_files has to be a list of fastq files (gzip supported).
    Output will be written to directory of first fastq file if outdir
    is not given. Assuming Sanger format of fastq files unless phred33
    is False.

    FIXME output file name (at least for paired): PDH132_s1.fastq.gz -> PDH132_s1_val_2.fq.gz

    FIXME Allow use of other trimmers
    """

    assert isinstance(fastq_files, list)

    try:
        trim_galore = CFG['trim_galore']
    except KeyError:
        trim_galore = None
    assert trim_galore

    cmd = [trim_galore]

    assert len(fastq_files) == 2# FIXME hardcoded paired
    cmd.extend(['--gzip', '--paired'])
    if not phred33:
        cmd.append('--phred64')
    if not outdir:
        outdir = os.path.dirname(fastq_files[0])
    cmd.extend(['-o', outdir])
    for f in fastq_files:
        cmd.append(f)

    # trim_galore also creates separate log
    #
    # FIXME missing simul opt
    # part of the problem is that the output files can't be specified
    (stdout, stderr) = run_cmd(cmd)

    # output files names can't specified. could be trimmed.fq if
    # single or val_X f paired and could be gzipped or not depending
    # on input and output settings. try to parse from output. if
    # --paired was given the files will be last element in lines
    # starting with 'Writing validated paired-end read'
    outfiles = []
    assert '--paired' in cmd
    for line in stderr:
        if line.startswith("Writing validated paired-end read"):
            f = line.split(' ')[-1]
            f = os.path.join(outdir, f)
            outfiles.append(f)
    return outfiles


def run_fastqc(fastq_files):
    """Run FastQC
    """

    assert isinstance(fastq_files, list)

    try:
        fastqc = CFG['fastqc']
    except KeyError:
        fastqc = None
    assert fastqc

    cmd = [fastqc]
    if CFG['java'] != 'java':
        cmd.extend(['-j', CFG['java']])
    cmd.extend(['-t', '%s' % CFG['numthreads']])
    cmd.extend(fastq_files)

    # FIXME missing simul opt.
    # needs to define output file name first.
    run_cmd(cmd)

# ----------------------------------------------------------------------
# EO unused bits
# ----------------------------------------------------------------------


def run_markdups(bam, outbam):
    """Run Picard's Markduplicate
    """

    markdups_jar = CFG['picard_dir'] + 'MarkDuplicates.jar'
    assert os.path.exists(markdups_jar)

    cmd = [CFG['java']]
    cmd.extend(CFG['java_opts'])
    cmd.extend(['-jar', markdups_jar])

    outmetrics = outbam.replace(".bam", "") + ".metrics"
    cmd.append("VALIDATION_STRINGENCY=LENIENT")
    cmd.append("INPUT=%s" % bam)
    cmd.append("TMP_DIR=%s" % os.path.dirname(outbam))
    cmd.append("OUTPUT=%s" % outbam)
    cmd.append("METRICS_FILE=%s" % outmetrics)
    cmd.append("CREATE_INDEX=true")

    if not CFG['simul']:
        log_fh = open("%s.log" % outbam, 'w')
        run_cmd(cmd, log_stdout=log_fh, log_stderr=log_fh)
        log_fh.close()
    else:
        assert not os.path.exists(outbam)# superflous? harmful? had an assert trigger if this was a link. removed assert and link was overwritten in simul mode
        fh = open(outbam, 'w')
        fh.write("cmd=%s\n" % ' '.join(cmd))
        fh.close()


def run_index_bam(bam):
    """Index BAM file with samtools

    FIXME untested
    """

    try:
        samtools = CFG['samtools']
    except KeyError:
        samtools = None
    assert samtools

    cmd = [samtools, 'index', bam]
    run_cmd(cmd)


def run_indelrealigner_targets(bam, ivals):
    """Run GATK's RealignerTargetCreator
    See http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_indels_RealignerTargetCreator.html
    """

    cmd = [CFG['java']]
    cmd.extend(CFG['java_opts'])
    cmd.extend(['-jar', CFG['gatk']])
    cmd.extend(['-nt', '%s' % CFG['numthreads']])
    cmd.extend(['-T', 'RealignerTargetCreator'])
    cmd.extend(['-I', bam])
    cmd.extend(['-R', CFG['reffa']])
    cmd.extend(['-o', ivals])
    if CFG['dbsnp']:
        cmd.extend(['--known', CFG['dbsnp']])

    if not CFG['simul']:
        log_fh = open("%s.log" % ivals, 'w')
        run_cmd(cmd, log_stdout=log_fh, log_stderr=log_fh)
        log_fh.close()
    else:
        assert not os.path.exists(ivals)
        fh = open(ivals, 'w')
        fh.write("cmd=%s\n" % ' '.join(cmd))
        fh.close()


def run_indelrealigner(bam, ivals, realn_bam):
    """Run GATK's IndelRealigner
    See http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_indels_IndelRealigner.html
    """

    cmd = [CFG['java']]
    cmd.extend(CFG['java_opts'])
    cmd.extend(['-jar', CFG['gatk']])
    cmd.extend(['-T', 'IndelRealigner'])
    cmd.extend(['-R', CFG['reffa']])
    cmd.extend(['-I', bam])
    cmd.extend(['-targetIntervals', ivals])
    cmd.extend(['-o', realn_bam])
    # no threading possible
    if CFG['dbsnp']:
        cmd.extend(['--knownAlleles', CFG['dbsnp']])

    if not CFG['simul']:
        log_fh = open("%s.log" % realn_bam, 'w')
        run_cmd(cmd, log_stdout=log_fh, log_stderr=log_fh)
        log_fh.close()
    else:
        assert not os.path.exists(realn_bam)
        fh = open(realn_bam, 'w')
        fh.write("cmd=%s\n" % ' '.join(cmd))
        fh.close()


def run_baserecalibrator_table(bam, recal_table):
    """Run GATK's BaseRecalibrator

    See http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_bqsr_BaseRecalibrator.html
    """

    cmd = [CFG['java']]
    cmd.extend(CFG['java_opts'])
    cmd.extend(['-jar', CFG['gatk']])
    cmd.extend(['-nct', '%s' % CFG['numthreads']])
    cmd.extend(['-T', 'BaseRecalibrator'])
    cmd.extend(['-R', CFG['reffa']])
    cmd.extend(['-I', bam])
    # mandatory
    assert CFG['dbsnp']
    cmd.extend(['-knownSites', CFG['dbsnp']])
    cmd.extend(['-o', recal_table])

    if not CFG['simul']:
        log_fh = open("%s.log" % recal_table, 'w')
        run_cmd(cmd, log_stdout=log_fh, log_stderr=log_fh)
        log_fh.close()
    else:
        assert not os.path.exists(recal_table)
        fh = open(recal_table, 'w')
        fh.write("cmd=%s\n" % ' '.join(cmd))
        fh.close()


def run_baserecalibrator(bam, recal_table, recal_bam):
    """Run GATK's PrintReads for BaseRecalibrator

    See http://gatkforums.broadinstitute.org/discussion/44/base-quality-score-recalibration-bqsr
    """

    cmd = [CFG['java']]
    cmd.extend(CFG['java_opts'])
    cmd.extend(['-jar', CFG['gatk']])
    cmd.extend(['-nct', '%s' % CFG['numthreads']])
    cmd.extend(['-T', 'PrintReads'])
    cmd.extend(['-R', CFG['reffa']])
    cmd.extend(['-I', bam])
    cmd.extend(['-BQSR', recal_table])
    cmd.extend(['-o', recal_bam])

    if not CFG['simul']:
        log_fh = open("%s.log" % recal_bam, 'w')
        run_cmd(cmd, log_stdout=log_fh, log_stderr=log_fh)
        log_fh.close()
    else:
        assert not os.path.exists(recal_bam)
        fh = open(recal_bam, 'w')
        fh.write("cmd=%s\n" % ' '.join(cmd))
        fh.close()



# ------------------------------------------------------------
# ruffus tasks
#
# FIXME make suffices variables
#
# ------------------------------------------------------------


# starting point
def generate_parameters():
    """Returning BAM input file
    """
    for f in CFG['bam_in']:
        y = (f, f.replace(".bam", ".mdups.bam"))
        yield y

@files(generate_parameters)
def markdups(input, output):
    """Proxy to run_markdups()
    """
    # FIXME shouldn't need to check but had a case where I pretended via links and touched files that mdup had already run but ruffus inisted on running it again!? posted to the ruffus group in 2014-03-26
    if not os.path.exists(output):
        run_markdups(input, output)


@follows(markdups)
# not arbitrary: GATK only recognizes certain extensions
@transform(markdups, suffix('.bam'), '.realn.intervals')
def indelrealigner_targets(input, output):
    """Proxy to run_indelrealigner_targets()
    """
    run_indelrealigner_targets(input, output)


@follows(indelrealigner_targets)
@transform(indelrealigner_targets, suffix('.realn.intervals'), add_inputs(r'\1.bam'), '.realn.bam')
def indelrealigner(input, output):
    """Proxy to run_indelrealigner()
    """
    run_indelrealigner(input[1], input[0], output)


@follows(indelrealigner)
@transform(indelrealigner, suffix('.realn.bam'), '.realn.recal.table')
def baserecalibrator_table(input, output):
    """Proxy to run_baserecalibrator_table()
    """
    run_baserecalibrator_table(input, output)


@follows(baserecalibrator_table)
@transform(baserecalibrator_table, suffix('.realn.recal.table'), add_inputs(r'\1.realn.bam'), '.realn.recal.bam')
def baserecalibrator(input, output):
    """Proxy to run_baserecalibrator_table()
    """
    run_baserecalibrator(input[1], input[0], output)


# ------------------------------------------------------------
# user interface
# ------------------------------------------------------------


def cmdline_parser():
    """
    creates an OptionParser instance
    """

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        dest="verbose",
                        help="be verbose")
    parser.add_argument("--debug",
                        action="store_true",
                        dest="debug",
                        help="enable debugging")
    parser.add_argument("--threads",
                        type=int,
                        dest="numthreads",
                        default=CFG['numthreads'],
                        help="Number of threads to use for calls"
                        " (default %d)" % CFG['numthreads'])
    default = 1
    parser.add_argument("--multiproc",
                        type=int,
                        default=default,
                        dest="nummultiproc",
                        help="Number of processes to run in parallel"
                        " (default %d)" % default)
    #parser.add_argument('-c', "--config",
    #                    required=True,
    #                    dest="config",
    #                    help="Config file")
    parser.add_argument('-i', "--bam",
                        required=True,
                        dest="bam",
                        nargs='+',
                        help="BAM input file/s")
    parser.add_argument('-r', "--reffa",
                        required=True,
                        dest="reffa",
                        help="Reference fasta file (indexed and with dict!"
                        " For this you can for example use Picard's CreateSequenceDictionary with CREATE_INDEX=true")
    parser.add_argument('-k', "--known",
                        required=True,
                        dest="dbsnp",
                        help="VCF file of known variants (for indel-realignment and base-call quality recalibration). Create an empty one just containing a header if not available")
                        # FIXME create empty one if arg is empty and warn (remove required True)
    parser.add_argument('-t', "--tasks",
                        dest="tasks",
                        default=DEFAULT_TASK,
                        help="Task to run (default %s)" % (DEFAULT_TASK))
    # FIXME output list of all tasks (currently not possible with ruffus)
    parser.add_argument("--only-print",
                        action="store_true",
                        dest="only_print",
                        help="Only print what would happen")
    parser.add_argument("--simul",
                        action="store_true",
                        dest="simul",
                        help="Only simulate execution. Creates fake output"
                        " files which only contain the command to be executed")
    return parser


def main():
    """main function
    """

    parser = cmdline_parser()
    args = parser.parse_args()

    ruffus_verbosity = 1
    if args.verbose:
        LOG.setLevel(logging.INFO)
        ruffus_verbosity = 2
    if args.debug:
        LOG.setLevel(logging.DEBUG)
        ruffus_verbosity = 4

    #with open(args.config) as configfile:
    #    config=json.load(configfile)
    CFG['bam_in'] = args.bam
    if args.numthreads:
        CFG['numthreads'] = args.numthreads
    # control number of threads in java's garbage collection
    # FIXME should be part of config
    CFG['java_opts'].append('-XX:ParallelGCThreads=%d' % CFG['numthreads'])
    if args.dbsnp:
        CFG['dbsnp'] = args.dbsnp

    CFG['reffa'] = args.reffa
    # FIXME the following two should be part of the pipeline
    faidx = CFG['reffa'] + ".fai"
    if not os.path.exists(faidx):
        LOG.fatal("faidx for reference fa %s missing. Please run samtools faidx first" % CFG['reffa'])
        sys.exit(1)
    fadict = os.path.splitext(CFG['reffa'])[0] + ".dict"
    if not os.path.exists(fadict):
        LOG.fatal("dict for reference fa %s missing. Please run Picard's CreateSequenceDictionary first" % CFG['reffa'])
        sys.exit(1)

    CFG['simul'] = args.simul


    #import pdb; pdb.set_trace()
    if args.only_print:
        pipeline_printout(sys.stdout, [args.task])
    else:
        pipeline_run([args.task], #touch_files_only=True,
                     multiprocess=args.nummultiproc,
                     verbose=ruffus_verbosity, logger=LOG)


if __name__ == "__main__":
    main()
    #LOG.info("Successful program exit")
    # FIXME add test set (E.coli or viral? Needs RG)






