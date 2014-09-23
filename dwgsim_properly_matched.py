#!/usr/bin/env python
"""Parse a BAM file containing reads simulated with DWGSIM and only
keep "correctly" mapped reads
"""


__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2014 Genome Institute of Singapore"
__license__ = "WTFPL http://www.wtfpl.net/"


import os
import sys

import pysam


BEEP_EVERY = 100000


def properly_matched(read_from, write_to):
    """FIXME
    """
    sam_in = pysam.Samfile(read_from, "rb" )
    sam_out = pysam.Samfile(write_to, "wb", template=sam_in)
    
    nwritten = 0
    for (nparsed, read) in enumerate(sam_in):
        if nparsed % BEEP_EVERY == 1:
            sys.stderr.write("Analyzing read %d\n" % nparsed)

        if not read.is_proper_pair:
            continue
    
        dwgsim = dict(zip(
            ("rid", "chrom", "startend1", "startend2", 
             "strandend1", "strandend2", "randend1", "randend2", 
             "stats1", "stats2", "rhash"),
            read.qname.split("_")))
        for k in ["startend1", "startend2", "strandend1", "strandend2", "randend1", "randend2"]:
            dwgsim[k] = int(dwgsim[k])
    
        chrom = sam_in.getrname(read.tid)
        if chrom == dwgsim["chrom"]:
            # lazy evaluation. can't figure out where (1 or 2?) pos is
            # stored and if we need to take care of complement (let alone
            # clipping)
            #
            # mate_no = dwgsim.rid[-1]
            write = False
            for k in ["startend1", "startend2"]:
                diff = abs(dwgsim[k]-read.pos)
                if diff < read.rlen:
                    write = True
                    break
            if write:
                nwritten += 1
                sam_out.write(read)
    sys.stderr.write("Wrote %d reads\n" % nwritten)
    sam_in.close()
    sam_out.close()
    

if __name__ == "__main__":
 
    try:
        read_from = sys.argv[1]
    except IndexError:
        read_from = "-"
    try:
        write_to = sys.argv[2]
    except IndexError:
        write_to = "-"

    if read_from != "-":
        assert os.path.exists(read_from), (
            "input file %s does bot exist" % read_from)
    if write_to != "-":
        assert not os.path.exists(write_to), (
            "refusing to overwrite already existing file %s" % write_to)
    sys.stderr.write("Reading from %s\n" % read_from)
    sys.stderr.write("Writing to %s\n" % write_to)

    properly_matched(read_from, write_to)
