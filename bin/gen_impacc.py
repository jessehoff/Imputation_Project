import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import seaborn as sb
%matplotlib inline
from sys import argv

script, imputed_gen, imputed_sample, true_gen, true_sample, impacc_log = argv

impfile = {}
with open(imputed_gen, 'r') as f:
    for line in f:
    #hapra=line.split(" ")[3:5]
        imppos=line.split(" ")[2]
        imp=[float(i) for i in (line.split(" ")[5:])]
        implist=[list(x) for x in zip(*[iter(imp)]*3)]
        altimp = [(x[1]+x[2]*2) for x in implist]
        impfile[imppos] = altimp
with open(imputed_sample, 'r') as f:
    impsamp = f.readlines()
impids = [item.split(" ")[1] for item in impsamp[2:]]
impdf = pd.DataFrame.from_dict(impfile, orient='index')
impdf.columns = impids

truefile = {}
with open(true_gen, 'r') as f:
    for line in f:
        refalt=line.split(" ")[3:5]
        pos=line.split(" ")[2]
        gen=[float(i) for i in (line.split(" ")[5:])]
        genlist=[list(x) for x in zip(*[iter(gen)]*3)]
        altgen = [(x[1]+x[2]*2) for x in genlist]
        truefile[pos] = altgen
with open(true_sample, 'r') as f:
    truesamp = f.readlines()
genids = [item.split(" ")[1] for item in truesamp[2:]]
truedf = pd.DataFrame.from_dict(truefile, orient='index')
truedf.columns = genids

corrs = impdf.corrwith(truedf, axis=1)
cor = pd.DataFrame(abs(corrs))
cor = cor.rename(index=str, columns={0:"correlation"})
mean_cor = str(cor.mean()).split("\n")[0].split("   ")[1]

with open(impacc_log, "a") as log:
    log.write(imputed_gen.split("/")[-1] + "\t" + mean_cor + "\n")
