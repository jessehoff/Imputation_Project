import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import seaborn as sb
import matplotlib
plt.style.use('fivethirtyeight')
from sys import argv

script, h, outfile = argv

assay= pd.read_table(h,delim_whitespace=True)
low_p = assay[(assay.P<=0.001)]
#Here, we're looking only at P-values below 10^-3
#This should make the histogram a little more interesting
ax = low_p.P.plot.hist(bins =50, color='r')
ax.set_axis_bgcolor('white')
ax.set_xlabel("HWE P-Value")
ax.set_ylabel("Number of SNPs")
z=h.strip("../hwe_stats/")
title=z.strip('.hwe')
ax.set_title(title + ' HWE P-Values')
plt.xscale('log')
for item in(ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(10)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
    item.set_fontsize(15)
fig = ax.get_figure()
fig.savefig(outfile, format='png')
plt.close('all')
