#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("-m", "--mpileup", help="input mpipleup, must be run with -r and --output-QNAME",  type=argparse.FileType('r'))
parser.add_argument("-n", "--nucfreq", help="input nucfreq file",  type=argparse.FileType('r') )
parser.add_argument("outfile",nargs="?", help="output matrix", type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument('-d', action="store_true", default=False)
parser.add_argument('--header', action="store_true", default=False)
args = parser.parse_args()

import pysam
import pandas as pd
import re

poses = []
refs = set()
for line in args.nucfreq:
	tokens = line.strip().split()
	poses.append(int(tokens[1]))	
	refs.add(tokens[0])

assert len(refs) == 1
ref = str(refs.pop())
#print(ref, file=sys.stderr)



def parseopt(opt, numseqs):
	matches=re.findall(".*(\+|\-)([0-9]+)([ACGTNacgtn]+).*", opt)
	replace = ""
	if len(matches) > 0:
		return( ["n"] * numseqs)
		# tried something more fancy but it does not work 
		for insdel, length, chars in matches:
			length = int(length)
			new = opt.replace("{}{}{}".format(insdel, length, chars[0:length]), "")
			print(opt, insdel, length, chars, new)
			opt = new
	else:
		out = []
		for char in opt:
			if char in [",", "."]:
				out.append( ".")
			else:
				out.append( "1")
		assert len(out) == numseqs, "{} {}".format(out, opt)
		return(out)



out = {}
allnames = set()

for line in args.mpileup:
	tokens = line.strip().split()
	conitg = tokens[0]
	pos = int(tokens[1]) 
	if(pos in poses):
		rbase = tokens[2].upper()
		numseqs = int(tokens[3])
		opts = parseopt(tokens[4], numseqs) 
		names = tokens[6].split(",")
		

		for i in range(numseqs):
			opt = opts[i]
			name = names[i] 
			allnames.add(name)
			if name not in out:
				out[name] = {}
			if pos not in out[name]:
				out[name][pos] = {}
			out[name][pos] = opt   

for name in allnames:
	cur = out[name]
	write = ""
	for pos in poses:
		if pos not in cur:
			write += "n"
		else:
			write += cur[pos]
	write += "\t" + name + "\n"
	args.outfile.write(write)




exit()







