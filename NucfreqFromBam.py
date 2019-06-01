#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("infile", nargs="?", help="input bam file",  type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument("outfile",nargs="?", help="output bam file", type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument('-d', action="store_true", default=False)
parser.add_argument('-a', help="output all positions", action="store_true", default=False)
parser.add_argument('--header', action="store_true", default=False)
args = parser.parse_args()

import numpy as np
import matplotlib.pyplot as plt
import pysam
import pandas as pd

bam = pysam.AlignmentFile(args.infile)
refs = {}
for read in bam.fetch(until_eof=True):
	ref = read.reference_name
	if(ref not in refs):
		if(args.a):
			refs[ref] = [0, 2147483648]
		else:
			refs[ref] = [2147483648, 0]
	
	start = read.reference_start
	end = read.reference_end
	#print(ref, start, end, file=sys.stderr)
	if(refs[ref][0] > start ):
		refs[ref][0]=start
	if(refs[ref][1] < end ):
		refs[ref][1]=end
	
#print(refs, file=sys.stderr)

def getCovByBase(contig, start, end):
	coverage = bam.count_coverage(contig, quality_threshold = None, start=start, stop=end)
	assert len(coverage) == 4
	cov = {}
	cov["A"] = coverage[0]
	cov["C"] = coverage[1]
	cov["T"] = coverage[2]
	cov["G"] = coverage[3]
	return(cov)

def makeTable(contig, cov):
	contiglen = len(cov["A"])
	out = []
	print("{} length {}".format(contig, contiglen) , file=sys.stderr)
	for i in range(0, contiglen):
		out.append( "{}\t{}\t{}\t{}\t{}\t{}\n".format(contig, i, cov["A"][i], cov["C"][i], cov["G"][i], cov["T"][i]))
	return(out)


#
# creates a table of nucleotide frequescies 
#
for contig in refs:
	start, end = refs[contig]
	print("getting coverage between {}:{}-{}".format(contig, start, end), file=sys.stderr)
	cov = getCovByBase(contig, start, end)
	contiglen = len(cov["A"])
	print("pysam coverage done for {}, used length {}".format(contig, contiglen), file=sys.stderr)
	if(contiglen > 0):
		for i in range(contiglen):
			args.outfile.write( "{}\t{}\t{}\t{}\t{}\t{}\n".format(contig, start + i, cov["A"][i], cov["C"][i], cov["G"][i], cov["T"][i]))
	
		
		#tbl = makeTable(contig, cov)
		#for line in tbl:
		#	args.outfile.write(line)







