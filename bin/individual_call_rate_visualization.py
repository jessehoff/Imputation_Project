#!/usr/bin/env python
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.style.use('fivethirtyeight')
from sys import argv

script, i, outfile = argv

infile= pd.read_table(i,delim_whitespace=True)

lines = infile.F_MISS
weights = np.ones_like(lines)/len(lines)
ax =infile.F_MISS.plot.hist(bins=50, ylim=(0,0.2),color='g', weights=weights)
ax.set_axis_bgcolor('white')
ax.set_xlabel("Proportion of Genotypes Missing per Individual")
ax.set_ylabel("Proportion of Individuals")
z=i.strip("../individual_stats/")
title=z.strip('.imiss')
ax.set_title(title + ' Individual Call Rates (n ='+ str(len(lines))+')')
for item in(ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(10)

for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
    item.set_fontsize(15)

fig = ax.get_figure()
fig.savefig(outfile,format='png')
plt.close('all')
