#!/usr/bin/env python
import pandas as pd
from sys import argv

script, hdfile,f250file, outfile = argv
print(hdfile,f250file,outfile)
hdmap = pd.read_table(hdfile, delim_whitespace=True,header=None,usecols=[2,3])
f250map =  pd.read_table(f250file, delim_whitespace=True,header=None,usecols=[2,3])
print(f250map.tail())
print(hdmap.tail())

joinmap = pd.concat([hdmap,f250map])

joinmap = joinmap.reindex(columns=[3,2])
joinmap = joinmap.sort_values([2])
joinmap=joinmap.drop_duplicates()
print(joinmap.tail())
joinmap.columns = ['position','cM']
joinmap = joinmap.sort_values(['position'])
joinmap['Genetic_Map(cM)'] = joinmap.position / 1000000
joinmap['COMBINED_rate(cM/Mb)'] = 1.0
outmap = joinmap.reindex(columns=['position','COMBINED_rate(cM/Mb)','Genetic_Map(cM)'])
outmap.iloc[-1,1] = 0

outmap.to_csv(outfile,sep=' ',index=False)

