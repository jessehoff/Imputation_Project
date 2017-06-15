#!/usr/bin/env python

import pandas as pd
import numpy as np
from itertools import islice
from sys import argv

script, chr, legend, sample1, sample2, newlegend, newsample = argv

legendlist = []
with open(legend, "r") as ref:		#This works on the concatenated legend of merged references file that is outputted by impute2
	head = list(islice(ref,1))
	lines = list(islice(ref,0, None))
	leg = [xx.strip().split() for xx in lines]
	for xx in leg:
		pra = "_".join(list(xx[1:4]))
		col = ":".join(list([chr,pra]))
		xx[0] = col
		xx = " ".join(xx)
		legendlist.append(xx)
legendlist.insert(0,head[0].strip())

with open(newlegend, "w") as rev:
    [rev.write(str(line) + "\n") for line in legendlist]



samples = []						#Works on sample files from eagle phasing (concatentates and reformats)
with open (sample1, 'r') as f:		#Needs to accept in order that they were given to impute2 -- will be in that order
	#head = list(islice(f, 0, 1))
	samp = list(islice(f, 2, None))	#Grabs only the sample lines
	for yy in samp:
		samples.append(yy.strip())
with open(sample2, 'r') as h:
	sampx = list(islice(h, 2, None))
	for ww in sampx:
		samples.append(ww.strip())
slist = [xx.split() for xx in samples]	#Appends to list which is then reordered and joined.  Written line-by-line to a new, revised sample file
writesamp = []
for entry in slist:
	xx = " ".join(list([entry[1], entry[0],entry[0], entry[2]]))
	writesamp.append(xx)
with open(newsample, "w") as new:
	[new.write(str(line) + "\n") for line in writesamp]
