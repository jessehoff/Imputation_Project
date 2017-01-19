#!/usr/bin/env python

import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import seaborn as sb
import matplotlib
plt.style.use('fivethirtyeight')
from sys import argv

script, a, outname = argv

assay= pd.read_table(a,delim_whitespace=True)
chromosomes = assay.NCHROBS.max() 

assay['Callrate']= assay['NCHROBS']/chromosomes

ax = assay.Callrate.plot.hist(bins = 100, ylim=(0,5000))
ax.set_axis_bgcolor('white')
ax.set_xlabel("Proportion of Individuals Missing SNP")
ax.set_ylabel("Number of SNPs")
z=a.strip("../allele_stats/")
title=z.strip('.frq')
ax.set_title(title + ' Allele Call Rates')
for item in(ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(10)

for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
    item.set_fontsize(15)

fig = ax.get_figure()
outfile = outname + '.individual_call_rate.png'
fig.savefig(outfile,format='png')
plt.close('all')
