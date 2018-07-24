#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("infile", nargs="?", help="input nucfreq file",  type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument("outdir",nargs="?", help="output direcotry for plots")
parser.add_argument('-d', action="store_true", default=False)
parser.add_argument('--header', action="store_true", default=False)
args = parser.parse_args()


import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import argparse 
import pandas as pd


df = pd.read_table(args.infile)
df.columns = ["contig", "position", "A", "C", "G", "T"]
sort = np.flip( np.sort(df[["A","C","G","T"]].values) , 1)
df["first"] = sort[:,0]
df["second"] = sort[:,1]
df["third"] = sort[:,2]
df["fourth"] = sort[:,3]

for contig, group in df.groupby(by="contig"):
	print(contig)
	
	truepos = group.position.values
	first = group["first"].values
	second = group["second"].values

	matplotlib.rcParams.update({'font.size': 18})
	fig, ax = plt.subplots( figsize=(16,9) )
	prime, = plt.plot(truepos, first, 'o', color="black", markeredgewidth=0.0, markersize=2, label = "most frequent base pair")
	sec, = plt.plot(truepos, second,'o', color="red",   markeredgewidth=0.0, markersize=2, label = "second most frequent base pair")
	ax.set_xlabel('Collapse Position (bp)')
	ax.set_ylabel('Sequence Read Depth')

	ylabels = [format(label, ',.0f') for label in ax.get_yticks()]
	xlabels = [format(label, ',.0f') for label in ax.get_xticks()]
	ax.set_yticklabels(ylabels)
	ax.set_xticklabels(xlabels)

	# Hide the right and top spines
	ax.spines["right"].set_visible(False)
	ax.spines["top"].set_visible(False)
	# Only show ticks on the left and bottom spines
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')

	lgnd = plt.legend(loc="upper right")
	for handle in lgnd.legendHandles:
		handle._sizes = ([300.0])


	plt.savefig(contig + ".png", dip=900)



