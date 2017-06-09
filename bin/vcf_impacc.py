import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import seaborn as sb
#%matplotlib inline
import glob as glob
from sys import argv
from itertools import islice

script, true, imputed, acc, out = argv
print(script,true,imputed,acc,out)

truedf = pd.read_pickle(true)
print(truedf.head())
genskey = {'./.':np.NaN, '0/0':0, '0/1':1 ,'1/1':2}

impdict = {}
with open(imputed, "r") as i:
	head = list(islice(i,5))
	impsamples = next(i).strip().split()[9:]
	for line in i:
		splat = line.strip().split()
		cp = (splat[0], splat[1])
		pos = ":".join(cp)
		gen = splat[9:]
		alt = [genskey[p] for p in gen]
		impdict[pos] = alt
impdf = pd.DataFrame.from_dict(impdict, orient='index')
print(impsamples[:10])
impdf.columns = [j.split('_')[1] for j in impsamples]

corrs=impdf.corrwith(truedf, axis = 1)
cor = pd.DataFrame(corrs)
cor = cor.rename(index=str, columns={0:"correlation"})
cor.to_csv(out)

# outfile = "./test_correlations.txt"
assay = imputed.split("/")[1].strip(".vcf")
correlation = cor.mean()[0]

with open (acc, 'a') as out:
	out.write(assay + "\t" + str(correlation)+"\n")
