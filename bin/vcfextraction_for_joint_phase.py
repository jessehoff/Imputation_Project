import pandas as pd
import random
import numpy as np
import glob
from collections import Counter
from collections import defaultdict
import sys

famfile = sys.argv[1]
run = famfile.split('.')[-3]


print(run,' run\n')
imputeanimal = []
fam = open(famfile)
for i in fam:
    imputeanimal.append(i.split()[1])
print(imputeanimal[:10])

samples = glob.glob('testset_assays/*list1.fam')
print(samples[1].split('/')[0].split('.')[0])
samples=sorted(samples, key=lambda sample: samples[1].split('/')[1].split('.')[0],reverse=True)


peranimal = defaultdict(list)
for idfile in samples: #open each output assay and figure out where those animals camefrom
    temp = []
    with open(idfile) as fp:
        for line in fp:
            name = line.split()[1]
            peranimal[name].append(idfile)

arrayset = defaultdict(list)
for i in imputeanimal:
    arrayset[peranimal[i][0]].append(i)

maps = glob.glob('./testset_assays/*.list1*.bim')
maps


for mapfile in maps: 
    name = 'merged_chrsplit/phased_' +  mapfile[17:-10] +'.' +run + '' + '.vcfregion'
    print(name,mapfile[17:-10])
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
    name = 'merged_chrsplit/phased_' +  k[15:-10]+'.list1.' +run + '' + '.keepvcf'
    print(name,k[15:-10])
    filout = open(name,'w')
    for i in v:
        filout.write(''.join([i,'\n']))
    filout.close()
