#!/usr/bin/env python
import pandas as pd
from sys import argv

script, infile, outfile = argv

buildmap = pd.read_table(infile, delim_whitespace=True,header=None,usecols=[2,3])
buildmap = buildmap.reindex(columns=[3,2])
buildmap.columns = ['position','cM']
buildmap['Genetic_Map(cM)'] = buildmap.position / 1000000
buildmap['COMBINED_rate(cM/Mb)'] = 1.0
outmap = buildmap.reindex(columns=['position','COMBINED_rate(cM/Mb)','Genetic_Map(cM)'])
outmap.iloc[-1,1] = 0

outmap.to_csv(outfile,sep=' ',index=False)

