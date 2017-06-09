#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('fivethirtyeight')

from sys import argv

script, corr, frq, map, hist, scatter, line, combo = argv

label = {'color':  'black',
		'weight': 'normal',
		'size': 10,
		}
title = {'color':  'black',
		'weight': 'normal',
		'size': 12,
		}

cor = pd.DataFrame.from_csv(corr)
mean = str(cor.mean()[0])

#Creates a simple histogram of correlations per base
hg = cor.correlation.plot.hist(bins = 100)
name = corr.strip(".snp_correlations.csv").strip("/imp_acc/")
hg.set_xlabel('True vs Imputed Correlation', fontdict=label)
hg.set_ylabel('Count', fontdict=label)
hg.grid(False)
hg.set_title("Histogram of Imputation Correlations per Locus\n" + name + '\nMean Correlation = ' + mean, fontdict=title)
fig = hg.get_figure()
fig.savefig(hist, format = 'png', bbox_inches = 'tight', transparent = True)
#plt.close('all')


#Data manipulation to merge MAF and correlation
freq = pd.read_table(frq, delim_whitespace=True)
del freq['A1']
del freq['A2']
del freq['NCHROBS']
map = pd.read_table(map, delim_whitespace=True, header=None)
map.columns = ("CHR", "SNP", "CM", "POS")
freq_map_merged = pd.merge(map, freq, on=("SNP","CHR"), how='inner')
freq_map_merged.head()
posmaf = freq_map_merged[["CHR","POS","MAF"]]
posmaf.index = posmaf.CHR.astype(str).str.cat(posmaf.POS.astype(str), sep = ':')
del posmaf.index.name
filter = posmaf["MAF"] <= 0.5
maffiltered = posmaf[filter]
mafcor = cor.join(maffiltered)

#
#Scatter plot of imputation MAF vs Imputation Accuracy (Correlation)
mafsp = mafcor.plot.scatter(x = "MAF", y = "correlation", alpha = 0.5)
name = corr.strip(".snp_correlations.csv").strip("/imp_acc/")
mafsp.set_xlabel('Minor Allele Frequency', fontdict=label)
mafsp.set_ylabel('True vs Imputed Correlation', fontdict=label)
mafsp.set_title("MAF vs Imputation Accuracy\n"+name, fontdict=title)
mafsp.grid(False)
fig = mafsp.get_figure()
fig.savefig(scatter, format = 'png', bbox_inches = 'tight', transparent = True)

#
#Putting Imp accuracies into bins
MAF = mafcor['MAF']
plotter = pd.DataFrame(MAF)
plotter["ACC"] = mafcor["correlation"]
bins = np.linspace(0,.5,100)
plotter_na = plotter[(-2<plotter['ACC']) & (plotter['ACC']<2)]
#print (plotter_na.describe())
a_bins = plotter_na.groupby(pd.cut(plotter_na['MAF'],bins))
mean_bin = a_bins.mean()
del mean_bin.index.name

#MAF vs Imputation Accuracy Line plot with binned values
mafline = a_bins.mean().plot(x='MAF',y='ACC', linewidth=2.0)
name = corr.strip(".snp_correlations.csv").strip("/imp_acc/")
mafline.set_xlabel('Minor Allele Frequency', fontdict=label)
mafline.set_ylabel('True vs Imputed Correlation', fontdict=label)
mafline.set_title("MAF vs Imputation Accuracy\n"+name,fontdict=title)
mafline.grid(False)
fig = mafline.get_figure()
fig.savefig(line, format = 'png', bbox_inches = 'tight', transparent = True)

plt.scatter(mafcor["MAF"], mafcor['correlation'], alpha = 0.1)
plt.plot(mean_bin['MAF'], mean_bin['ACC'])
plt.axis([0,0.5,0,1])
plt.grid(alpha = 0.0)
plt.ylabel('True vs. Imputed Genotype Correlation', fontsize=10, color='k')
plt.xlabel('Minor Allele Frequency', fontsize=10, color='k')
plt.title('Chip-level Imputation Accuracy by MAF\n'+ name + '\nMean Correlation = ' + mean, fontsize = 13)
plt.tick_params(colors='k')
plt.legend('')
plt.savefig(combo, format = 'png', bbox_inches = 'tight', transparent = True)

plt.close('all')
