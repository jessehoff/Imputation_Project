#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('fivethirtyeight')

from sys import argv
script, a, outfile = argv

assay= pd.read_table(a,delim_whitespace=True)
chromosomes = assay.NCHROBS.max() 

assay['Callrate']= assay['NCHROBS']/chromosomes

lines = assay.Callrate
weights = np.ones_like(lines)/len(lines)
ax = assay.Callrate.plot.hist(bins = 100, ylim=(0,0.1), weights=weights)
ax.set_axis_bgcolor('white')
ax.set_xlabel("Proportion of Individuals with SNP")
ax.set_ylabel("Number of SNPs")
plt.xticks(np.arange(min(assay.Callrate), max(assay.Callrate), .1))
z=a.strip("../snp_stats/")
title=z.strip('.frq')
ax.set_title(title + ' SNP Call Rates (n ='+ str(len(lines))+')') 
for item in(ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(10)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
    item.set_fontsize(12)

fig = ax.get_figure()
fig.savefig(outfile,format='png')
plt.close('all')
