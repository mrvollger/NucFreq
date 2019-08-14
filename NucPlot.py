#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("infile", nargs="?", help="input bam file",  type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument("outfile",nargs="?", help="output plot file")
parser.add_argument('-d', action="store_true", default=False)
parser.add_argument('-a', help="output all positions", action="store_true", default=False)
parser.add_argument('-r', '--repeatmasker', help="rm out to add to plot", type=argparse.FileType('r') , default=None)
parser.add_argument('-y', '--ylim', help="max y axis limit", type=float , default=None)
parser.add_argument('--header', action="store_true", default=False)
parser.add_argument("--psvsites", help="CC/mi.gml.sites", default=None)
parser.add_argument('-s', '--soft', action="store_true", default=False)
parser.add_argument('-c', '--minclip', help="min number of clippsed bases in order to be displayed", type=float , default=1000)
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



def getSoft(read):
	rtn = []
	cigar = read.cigartuples
	start = cigar[0]
	end = cigar[-1]
	if(start[0]==S):
		rtn.append( (read.reference_name, "start", start[1], read.reference_start, read) )
	if(end[0] == S):
		rtn.append( (read.reference_name, "end", end[1], read.reference_end, read) )
	return(rtn)


soft = []
bam = pysam.AlignmentFile(args.infile)
refs = {}
for read in bam.fetch(until_eof=True):
	ref = read.reference_name
	soft += getSoft(read)
	if(ref not in refs):
		if(args.a):
			refs[ref] = [0, 2147483648]
		else:
			refs[ref] = [2147483648, 0]
	
	start = read.reference_start
	end = read.reference_end
	if(refs[ref][0] > start ):
		refs[ref][0]=start
	if(refs[ref][1] < end ):
		refs[ref][1]=end
	


def getCovByBase(contig, start, end):
	coverage = bam.count_coverage(contig, quality_threshold = None, start=start, stop=end, read_callback="nofilter")
	assert len(coverage) == 4
	cov = {}
	cov["A"] = coverage[0]
	cov["C"] = coverage[1]
	cov["T"] = coverage[2]
	cov["G"] = coverage[3]
	return(cov)


#
# creates a table of nucleotide frequescies 
#
nf = []
for contig in refs:
	start, end = refs[contig]
	cov = getCovByBase(contig, start, end)
	contiglen = len(cov["A"])
	if(contiglen > 0):
		for i in range(contiglen):
			nf.append( [contig, start + i, cov["A"][i], cov["C"][i], cov["G"][i], cov["T"][i] ]   )
	
df = pd.DataFrame(nf, columns=["contig", "position", "A", "C", "G", "T"])
sort = np.flip( np.sort(df[["A","C","G","T"]].values) , 1)
df["first"] = sort[:,0]
df["second"] = sort[:,1]
df["third"] = sort[:,2]
df["fourth"] = sort[:,3]
df.sort_values(by=["contig", "position", "second"], inplace=True)


soft = pd.DataFrame(soft, columns=["contig", "side", "value", "position", "read"])
#print(len(soft), file=sys.stderr)
soft = soft[soft.value >= args.minclip]
#print(len(soft), file=sys.stderr)
soft.sort_values(by=["contig", "position"], inplace=True)



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



counter = 0
for contig, group in df.groupby(by="contig"):
	print(contig)
	
	truepos = group.position.values
	first = group["first"].values
	second = group["second"].values

	matplotlib.rcParams.update({'font.size': 18})
	fig, ax = plt.subplots( figsize=(16,9) )



	if(RM is not None):
		rmax = ax
		rm = RM[ (RM.qname == contig) & (RM.start >= min(truepos)) & (RM.end <= max(truepos)) ]
		assert len(rm.index) != 0, "No matching RM contig"
		#rmax.set_xlim(rm.start.min(), rm.end.max())
		#rmax.set_ylim(-1, 9)

		for idx, row in rm.iterrows():	
			width = row.end - row.start
			rect = patches.Rectangle((row.start,0), width, -max(first)/20, linewidth=1, edgecolor='none',facecolor=row.color, alpha = .75) 
			rmax.add_patch(rect)
		plt.show()

	prime, = ax.plot(truepos, first, 'o', color="black", markeredgewidth=0.0, markersize=2, label = "most frequent base pair")
	sec, = ax.plot(truepos, second,'o', color="red",   markeredgewidth=0.0, markersize=2, label = "second most frequent base pair")
	

	#inter = int( (max(truepos)-min(truepos))/50)
	#sns.lineplot(  (truepos/inter).astype(int)*inter, first, ax = ax, err_style="bars") 


	maxval = max(truepos)
	minval = max(truepos)
	if( maxval < 1000000 ):
		xlabels = [format(label, ',.0f') for label in ax.get_xticks()]
		lab = "bp"
	elif( maxval < 10000000):
		xlabels = [format(label/1000, ',.1f') for label in ax.get_xticks()]
		lab = "kbp"
	else:
		xlabels = [format(label/1000000, ',.2f') for label in ax.get_xticks()]
		lab = "Mbp"

	if(args.ylim is not None):
		ax.set_ylim(0, args.ylim)

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

	if(args.soft ):
		tmpsoft = soft[ soft.contig == contig ]
		if(len(tmpsoft) > 0 ):
			axsoft = ax.twinx()
			axsoft.invert_yaxis()
			bins = 150 
			color = "darkgreen"
			sns.distplot(tmpsoft.position, bins=bins, kde=False, ax=axsoft, hist_kws={'weights': tmpsoft.value/1000, "alpha": .25, "color": color})
			bot, top = axsoft.get_ylim()
			axsoft.set_ylim(1.1*bot, 0)
			# Hide the right and top spines
			axsoft.spines["top"].set_visible(False)
			# Only show ticks on the left and bottom spines
			axsoft.yaxis.set_ticks_position('right')
			axsoft.xaxis.set_ticks_position('bottom')
			axsoft.tick_params(axis='y', colors=color)
			axsoft.set_ylabel('Clipped Bases (kbp)', color=color)


	if(args.psvsites is not None):
		cuts = {}
		for idx, line in enumerate(open(args.psvsites).readlines()):
			try:	
				vals = line.strip().split()
				cuts[idx] = list(map(int, vals))
				#make plot
				x = np.array(cuts[idx]) -1
				idxs = (np.isin(truepos, x))
				y = second[ idxs ]
				ax.plot(x,y, alpha=0.5) #, label="group:{}".format(idx) )
			except Exception as e:
				print("Skipping because error: {}".format(e), file=sys.stderr)
				continue

	outpath = os.path.abspath(args.outfile)
	mydir = os.path.dirname(outpath) 
	name = os.path.basename(outpath)
	if(counter == 0):
		outf = 	outpath
	else:
		outf = "{}/{}_{}".format(mydir, counter, name)

	plt.savefig(outf, dip=900)
	counter += 1





