import pandas as pd
import numpy as np
import glob as glob
from itertools import islice
from sys import argv


script, true, truepickle= argv
print(script,true,truepickle)


truedict = {}
genskey = {'.|.':np.NaN, '0|0':0, '0|1':1 ,'1|1':2}
with open(true, "r") as t:
	head = list(islice(t,15))
	#print(head)
	samples = next(t).strip().split()[9:]
	for line in t:
		splat = line.strip().split()
		cp = (splat[0], splat[1])
		pos = ":".join(cp)
		gen = splat[9:]
		alt = [genskey[i] for i in gen]
		truedict[pos] = alt
truedf = pd.DataFrame.from_dict(truedict, orient='index')
truedf.columns = samples#[j.split('\t') for j in samples]
truedf.to_pickle(truepickle)
print(truedf.head())
