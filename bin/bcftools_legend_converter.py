#!/usr/bin/env python

import pandas as pd
import numpy as np
from itertools import islice
from sys import argv

script, chr, legend, newlegend = argv

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
