#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("bam", nargs="?", help="input bam file",  type=argparse.FileType('r'))
parser.add_argument("ref", nargs="?", help="input bam file",  type=argparse.FileType('r'))
parser.add_argument("outfile",nargs="?", help="output plot file")
parser.add_argument("-t", "--threshold", help="output plot file", type=int, default=80)
args = parser.parse_args()

import os 
import numpy as np
import pysam
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns 
from Bio import SeqIO
M=0 #M  BAM_CMATCH      0
I=1 #I  BAM_CINS        1
D=2 #D  BAM_CDEL        2
N=3 #N  BAM_CREF_SKIP   3
S=4 #S  BAM_CSOFT_CLIP  4
H=5 #H  BAM_CHARD_CLIP  5
P=6 #P  BAM_CPAD        6
E=7 #=  BAM_CEQUAL      7
X=8 #X  BAM_CDIFF       8
B=9 #B  BAM_CBACK       9
NM=10 #NM       NM tag  10
conRef  =       [M, D, N, E, E] # these ones "consume" the reference
conQuery=       [M, I, S, E, X] # these ones "consume" the query
conAln  =       [M, I, D, N, S, E, X] # these ones "consume" the alignments


bam = pysam.AlignmentFile(args.bam , "rb")
recs = list(SeqIO.parse(args.ref, "fasta"))

for rec in recs:
	ref = rec.name
	coverage = bam.count_coverage(ref, quality_threshold = None)
	coverage = {"A":coverage[0], "C":coverage[1], "G":coverage[2], "T":coverage[3] }
	df = pd.DataFrame.from_dict(coverage)
	df["contig"] = ref
	df["position"] = list(range(len(df)))

	#df = pd.DataFrame(nf, columns=["contig", "position", "A", "C", "G", "T"])
	sort = np.flip( np.sort(df[["A","C","G","T"]].values) , 1)
	df["first"] = sort[:,0]
	df["second"] = sort[:,1]
	df["third"] = sort[:,2]
	df["fourth"] = sort[:,3]
	df["cov"] = df["A"] + df["C"] + df["G"] + df["T"]

	print(df)
	assert len(df) == len(rec.seq)
	
	



