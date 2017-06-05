import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import seaborn as sb
#%matplotlib inline
import glob as glob
from sys import argv
from itertools import islice

script, true, imputed, log, out = argv

truedict = {}
genskey = {'./.':'NA', '0/0':'0', '0/1':'1' ,'1/1':'2'}
with open(true, "r") as t:
	head = list(islice(t,6))
	samples = t.next().strip().split()[9:]
	for line in t:
		splat = line.strip().split()
		cp = (splat[0], splat[1])
		pos = ":".join(cp)
		gen = splat[9:]
		alt = [genskey[i] for i in gen]
		truedict[pos] = alt
truedf = pd.DataFrame.from_dict(truedict, orient='index')
truedf.columns = samples

impdict = {}
with open(imputed, "r") as i:
	head = list(islice(i,6))
	impsamples = i.next().strip().split()[9:]
	for line in i:
		splat = line.strip().split()
		cp = (splat[0], splat[1])
		pos = ":".join(cp)
		gen = splat[9:]
		alt = [genskey[p] for p in gen]
		impdict[pos] = alt
impdf = pd.DataFrame.from_dict(impdict, orient='index')
impdf.columns = impsamples

corrs=impdf.corrwith(truedf, axis = 1)
cor = pd.DataFrame(corrs)
cor = cor.rename(index=str, columns={0:"correlation"})

cor.to_csv(out)

outfile = "./test_correlations.txt"
assay = imputed.split("/")[1].strip(".vcf")
correlation = cor.mean()[0]

with open (log, 'a') as out:
	out.write(assay + "\t" + str(correlation)+"\n")

accuracy_test/merged_refs/F250_HD_merged.chr25.vcf
