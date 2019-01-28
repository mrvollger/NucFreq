#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("infile", nargs="?", help="input nucfreq file",  type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument("outdir",nargs="?", help="output direcotry for plots", default=".")
parser.add_argument('-d', action="store_true", default=False)
parser.add_argument('--header', action="store_true", default=False)
parser.add_argument('--second', default=None)
args = parser.parse_args()


import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import argparse 
import pandas as pd
import seaborn as sns


df = pd.read_table(args.infile)
df.columns = ["contig", "position", "A", "C", "G", "T"]
sort = np.flip( np.sort(df[["A","C","G","T"]].values) , 1)
df["first"] = sort[:,0]
df["second"] = sort[:,1]
df["third"] = sort[:,2]
df["fourth"] = sort[:,3]

if(args.second is not None):
	df2 = pd.read_table(args.second)
	df2.columns = ["contig", "position", "A", "C", "G", "T"]
	sort = np.flip( np.sort(df2[["A","C","G","T"]].values) , 1)
	df2["first"] = sort[:,0]
	df2["second"] = sort[:,1]
	df2["third"] = sort[:,2]
	df2["fourth"] = sort[:,3]
	print( (df["second"]).describe() )
	print( (df2["second"]).describe() )
	
	df["ratio"] = df["second"]/(df["second"] + df2["second"])
	print(df["ratio"].describe())
	fig, ax = plt.subplots( figsize=(16,9) )
	tmp = df["ratio"].dropna()
	#plt.density(tmp, bins=25)
	#plt.savefig("hist.png", dip=900)
	x = sns.distplot(tmp, hist=True).get_figure()
	x.savefig("hist.png")
	plt.clf()
	sns.distplot(df["second"]+df2["second"], hist=True).get_figure()
	sns.distplot(df["second"], hist=True).get_figure()
	x = sns.distplot(df2["second"], hist=True).get_figure()
	x.savefig("hist2.png")

for contig, group in df.groupby(by="contig"):
	print(contig)
	
	truepos = group.position.values
	first = group["first"].values
	second = group["second"].values

	matplotlib.rcParams.update({'font.size': 18})
	fig, ax = plt.subplots( figsize=(16,9) )
	prime, = plt.plot(truepos, first, 'o', color="black", markeredgewidth=0.0, markersize=2, label = "most frequent base pair")
	sec, = plt.plot(truepos, second,'o', color="red",   markeredgewidth=0.0, markersize=2, label = "second most frequent base pair")
	if(args.second is not None):
		temp = df2["second"].values
		tri, = plt.plot(truepos, temp,'o', color="green",   markeredgewidth=0.0, 
				markersize=2, label = "second most frequent base pair (forward)")

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


	plt.savefig(args.outdir + "/" + contig + ".png", dip=900)



