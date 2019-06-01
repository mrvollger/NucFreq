#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("infile", nargs="?", help="input nucfreq file",  type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument("outdir",nargs="?", help="output direcotry for plots", default=".")
parser.add_argument('-d', action="store_true", default=False)
parser.add_argument('--header', action="store_true", default=False)
parser.add_argument('-r', '--repeatmasker', help="rm out to add to plot", type=argparse.FileType('r') , default=None)
parser.add_argument('--second', default=None)
args = parser.parse_args()


import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import argparse 
import pandas as pd
import seaborn as sns

RM=None
colors = sns.color_palette()
cmap = {}
counter = 0
if(args.repeatmasker is not None):
	names = ["score", "perdiv", "perdel", "perins", "qname", "start", "end", "left", "strand", "repeat", "family", "rstart", "rend", "rleft", "ID"]
	lines = []
	for idx, line in enumerate(args.repeatmasker):
		if(idx > 2):
			lines.append(line.strip().split()[0:15])
	
	RM = pd.DataFrame(lines, columns=names)
	RM.start = RM.start.astype(int)
	RM.end = RM.end.astype(int)
	RM["label"] =RM.family.str.replace("/.*", "")
	for idx, lab in enumerate(sorted(RM.label.unique())):
		cmap[lab] = colors[ counter%len(colors) ]
		counter += 1
	RM["color"] = RM.label.map(cmap)
	
	args.repeatmasker.close()


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
	
	if(RM is not None):
		fig, (rmax, ax) = plt.subplots(2,1, figsize=(16,9), gridspec_kw = {'height_ratios':[1, 15]})
		rm = RM[RM.qname == contig]
		rmax.set_xlim(rm.start.min(), rm.end.max())
		rmax.set_ylim(0, 1)
		rmax.tick_params(
				axis='both',	# changes apply to the x-axis
				which='both',   # both major and minor ticks are affected
				left =False, right=False, labelleft=False, 
				bottom=False,   # ticks along the bottom edge are off
				top=False,       # ticks along the top edge are off
				labelbottom=False)

		for idx, row in rm.iterrows():	
			width = row.end - row.start
			rect = patches.Rectangle((row.start,0), width, 1, linewidth=1, edgecolor='none',facecolor=row.color, alpha = .75) 
			rmax.add_patch(rect)
		plt.show()
	else:
		fig, ax = plt.subplots( figsize=(16,9) )

	prime, = ax.plot(truepos, first, 'o', color="black", markeredgewidth=0.0, markersize=2, label = "most frequent base pair")
	sec, = ax.plot(truepos, second,'o', color="red",   markeredgewidth=0.0, markersize=2, label = "second most frequent base pair")
	if(args.second is not None):
		temp = df2["second"].values
		tri, = plt.plot(truepos, temp,'o', color="green",   markeredgewidth=0.0, 
				markersize=2, label = "second most frequent base pair (forward)")
	
	maxval = max(truepos)
	minval = max(truepos)
	if( maxval < 1000000 ):
		xlabels = [format(label, ',.0f') for label in ax.get_xticks()]
		lab = "bp"
	elif( maxval < 10000000):
		xlabels = [format(label/1000, ',.0f') for label in ax.get_xticks()]
		lab = "kbp"
	else:
		xlabels = [format(label/1000000, ',.1f') for label in ax.get_xticks()]
		lab = "Mbp"


	ax.set_xlabel('Collapse Position ({})'.format(lab))
	ax.set_ylabel('Sequence Read Depth')

	ylabels = [format(label, ',.0f') for label in ax.get_yticks()]
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



