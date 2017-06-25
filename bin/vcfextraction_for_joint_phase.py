import pandas as pd
import random
import numpy as np
import glob
from collections import Counter
from collections import defaultdict
import sys
import re

famfile = sys.argv[1]
nlist = sys.argv[2]
run = famfile.split('/')[1][3:]


print(run,' run\n')
imputeanimal = []
fam = open(famfile)
for i in fam:
    imputeanimal.append(i.split()[1])
print(imputeanimal[:10])

samples = glob.glob('downsample/*list'+nlist + '.fam')
print(samples)
print(samples[1].split('/')[0].split('.')[0])
samples=sorted(samples, key=lambda sample: samples[1].split('/')[1].split('.')[0],reverse=True)


peranimal = defaultdict(list)
for idfile in samples: #open each output assay and figure out where those animals camefrom
    temp = []
    with open(idfile) as fp:
        for line in fp:
            name = line.split()[1] #this is a 
            peranimal[name].append(idfile)


arrayset = defaultdict(list)
for i in imputeanimal:
    arrayset[peranimal[i][0]].append(i) #only appends the first time the animal appears

maps = glob.glob('./downsample/*.list'+nlist+'*.bim')

mapassay = re.search(r'(\ggpld).list', maps[1])

for mapfile in maps: 
    match = re.search(r'(\w*).list', mapfile)
    mapname = match.group(1)
    name = 'merged_chrsplit/run'+run+'/phased_' +  mapname   + '.vcfregion'
    filout = open(name,'w')
    with open(mapfile) as fp:
        for line in fp:
            site = line.split()[3]
            chr = line.split()[0]
            filout.write(chr + '\t' + site + '\n')
    filout.close()


#merged_chrsplit/phased_snp50.list1.run1.keepvcf
#merged_chrsplit/phased_snp50.list1.run1.vcfregion


for k,v in arrayset.items():
    match = re.search(r'(\w*).list', k)
    mapname = match.group(1)
    name = 'merged_chrsplit/run'+run+'/phased_' +  mapname+'.keepvcf'
    print(name,mapname)
    filout = open(name,'w')
    for i in v:
        filout.write(''.join([i,'\n']))
    filout.close()
