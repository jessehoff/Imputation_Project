#!/usr/bin/env python
import datetime
import time
import pandas as pd
import numpy as np
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sb
from sys import argv

script, corrpath, frq, mapfile, corrout, fig, lowmaf, masterstats = argv

name = str(corrpath.split('/')[-1].split('.')[0])
run = str(corrpath.split('/')[1])

interesting_files = glob.glob(corrpath+"*.csv")
df_list = []
for filename in sorted(interesting_files):
    df_list.append(pd.read_csv(filename))
full_df = pd.concat(df_list)
cor = full_df
cor.index = cor['Unnamed: 0']
del cor['Unnamed: 0']
del cor.index.name
cor['correlation']  = abs(cor['correlation'])

total_sites = str(len(cor))

freq = pd.read_table(frq, delim_whitespace=True)
del freq['A1'] #you can use use_cols to pick which columsn to read in to start. 
del freq['A2']
del freq['NCHROBS']
map = pd.read_table(mapfile, delim_whitespace=True, header=None)
map.columns = ("CHR", "SNP", "CM", "POS")
freq_map_merged = pd.merge(map, freq, on=("SNP","CHR"), how='inner')
freq_map_merged.head()
posmaf = freq_map_merged[["CHR","POS","MAF"]]
posmaf.index = posmaf.CHR.astype(str).str.cat(posmaf.POS.astype(str), sep = ':')
del posmaf.index.name
filter = posmaf["MAF"] <= 0.5
maffiltered = posmaf[filter]
mafcor = cor.join(maffiltered)

MAF = mafcor['MAF']
plotter = pd.DataFrame(MAF)
plotter["ACC"] = mafcor["correlation"]
bins = np.linspace(0,.5,100)
plotter_na = plotter[(-2<plotter['ACC']) & (plotter['ACC']<2)]
a_bins = plotter_na.groupby(pd.cut(plotter_na['MAF'],bins))
mean_bin=a_bins.mean()
del mean_bin.index.name

impmean = str(mafcor['correlation'].mean())

fixed = mafcor[mafcor.correlation.isnull()&mafcor.MAF.notnull()]
nfixed = str(fixed['MAF'].count())
fixedmaf = str(fixed['MAF'].mean())


corrout = mafcor.to_csv(corrout)

plt.style.use('ggplot')
plt.scatter(mafcor["MAF"], mafcor['correlation'], alpha = 0.05)
plt.plot(mean_bin['MAF'], mean_bin['ACC'], color = 'r')
plt.axis([0,0.5,0,1])
plt.grid(alpha = 0.0)
plt.ylabel('True vs. Imputed Genotype Correlation', fontsize=10, color='k')
plt.xlabel('Minor Allele Frequency', fontsize=10, color='k')
plt.title('900K Imputation Accuracy by MAF\n' + name + ' Imputed - '+ run)
plt.tick_params(colors='k')
plt.savefig(fig)

#Performing same analysis with MAF < 0.10

MAF = mafcor['MAF']
plotter = pd.DataFrame(MAF)
plotter["ACC"] = mafcor["correlation"]
bins = np.linspace(0,.1,100)
plotter_na = plotter[(-2<plotter['ACC']) & (plotter['ACC']<2)]
a_bins = plotter_na.groupby(pd.cut(plotter_na['MAF'],bins))
filter = mafcor['MAF'] <= 0.1
mafcor = mafcor[filter]

nlowmaf = str(mafcor['MAF'].count())
lowmafimpacc = str(mafcor['correlation'].mean())

plt.style.use('ggplot')
plt.scatter(mafcor["MAF"], mafcor['correlation'], alpha = 0.05)
plt.plot(mean_bin['MAF'], mean_bin['ACC'], color = 'r')
plt.axis([0,0.1,0,1])
plt.grid(alpha = 0.0)
plt.ylabel('True vs. Imputed Genotype Correlation', fontsize=10, color='k')
plt.xlabel('Minor Allele Frequency', fontsize=10, color='k')
plt.title('900K Imputation Accuracy by MAF (< 0.10)\n'+ name + ' Imputed - '+ run)
plt.tick_params(colors='k')
plt.savefig(lowmaf)

ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

x = [run, name, total_sites, impmean, nfixed, fixedmaf, nlowmaf, lowmafimpacc, st]
line = '\t'.join(x)

with open(masterstats, 'a') as m:
	m.write(line + '\n')
