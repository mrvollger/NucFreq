#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("infile", nargs="?", help="input bam file",  type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument("outfile",nargs="?", help="output bam file", type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument('-d', action="store_true", default=False)
parser.add_argument('--header', action="store_true", default=False)
args = parser.parse_args()

import numpy as np
import matplotlib.pyplot as plt
import pysam
import pandas as pd

bam = pysam.AlignmentFile(args.infile)
refs = bam.references


def getCovByBase(contig):
	coverage = bam.count_coverage(contig, quality_threshold = None)
	assert len(coverage) == 4
	cov = {}
	cov["A"] = coverage[0]
	cov["C"] = coverage[1]
	cov["T"] = coverage[2]
	cov["G"] = coverage[3]
	return(cov)

def makeTable(contig, cov):
	contiglen = len(cov["A"])
	out = ""
	for i in range(0, contiglen):
		out += "{}\t{}\t{}\t{}\t{}\t{}\n".format(contig, i, cov["A"][i], cov["C"][i], cov["G"][i], cov["T"][i])
	return(out)


#
# creates a table of nucleotide frequescies 
#
nucfreq = ""
for contig in refs:
	cov = getCovByBase(contig)
	tbl = makeTable(contig, cov)
	nucfreq += tbl
args.outfile.write(nucfreq)






