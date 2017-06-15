#!/usr/bin/env python

import pandas as pd
import numpy as np
from itertools import islice
from sys import argv

script, true, imputed, acc, out = argv

truedf = pd.read_pickle(true)
impdict = {}

with open(imputed, 'r') as i:
    head = list(islice(i,10))
    impsamples = next(i).strip().split()[9:]
    for line in i:
        splat = line.strip().split()
        pos = splat[2]
        gen = splat[9:]
        alt = [float(xx.split(":")[1]) for xx in gen]
        impdict[pos] = alt
impdf = pd.DataFrame.from_dict(impdict, orient='index')
impdf.columns = impsamples

corrs=impdf.corrwith(truedf, axis = 1)
cor = pd.DataFrame(corrs)
cor = cor.rename(index=str, columns={0:"correlation"})
cor.to_csv(out)

assay = imputed.split("/")[1].strip(".dose.vcf")
correlation = cor.mean()[0]
with open (acc, 'a') as out:
	out.write(assay + "\t" + str(correlation)+"\n")
