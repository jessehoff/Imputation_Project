import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import seaborn as sb
#%matplotlib inline
import glob as glob
from sys import argv

script, true, imputed, log = argv

truedict = {}
trueids = []
with open(true, "r") as t:
    for line in t:
        if not line.startswith("#"):
            cp = (line.split("\t")[0], line.split("\t")[1])       #
            pos = ":".join(cp)
            gen = line.split("\t")[9:]
            genotypes = [(x.split("/")) for x in gen]
            intgt = [[int(y) for y in x if y != '.' and y != '.\n']for x in genotypes]
            alt = [sum(x) for x in intgt]
            truedict[pos] = alt
        if line.startswith("#"):
            trueids.append(line)
samples = trueids[-1].split("\t")[9:]
truedf = pd.DataFrame.from_dict(truedict, orient='index')
truedf.columns = samples

impdict = {}
impids = []
with open(imputed, "r") as i:
    for line in i:
        if not line.startswith("#"):
            cp = (line.split("\t")[0], line.split("\t")[1])       #
            pos = ":".join(cp)
            gen = line.split("\t")[9:]
            genotypes = [(x.split("/")) for x in gen]
            intgt = [[int(y) for y in x if y != '.' and y != '.\n']for x in genotypes]
            alt = [sum(x) for x in intgt]
            impdict[pos] = alt
        if line.startswith("#"):
            impids.append(line)
impsamples = impids[-1].split("\t")[9:]
impdf = pd.DataFrame.from_dict(impdict, orient='index')
impdf.columns = impsamples

corrs=impdf.corrwith(truedf, axis = 1)
cor = pd.DataFrame(corrs)
cor = cor.rename(index=str, columns={0:"correlation"})

outfile = "./test_correlations.txt"
assay = imputed.split("/")[1].strip(".vcf")
correlation = cor.mean()[0]

with open (log, 'a') as out:
    out.write(assay + "\t" + str(correlation)+"\n")
