import pandas as pd
import random
import numpy as np
import glob
from collections import Counter
from collections import defaultdict
import sys

famfile = sys.argv[1]

imputeanimal = []
fam = open(famfile)
for i in fam:
    imputeanimal.append(i.split()[1])
print(imputeanimal[:10])

samples = glob.glob('../dataprepper/testset_data/*hol_testset*.fam')
print(samples[1].split('/')[3].split('.')[2])
samples=sorted(samples, key=lambda sample: int(sample.split('/')[3].split('.')[2]),reverse=True)
print(samples[:10])

peranimal = defaultdict(list)
for idfile in samples:
    temp = []
    with open(idfile) as fp:
        for line in fp:
            name = line.split()[1]
            peranimal[name].append(idfile)

arrayset = defaultdict(list)
for i in imputeanimal:
    arrayset[peranimal[i][0]].append(i)

maps = glob.glob('../dataprepper/testset_data/*hol_testset*.bim')
maps


for mapfile in maps:
    name = 'phased_' +  mapfile.lstrip("../dataprepper/testset_data/").rstrip('.bim') + '.vcfregion'
    print(name)
    filout = open(name,'w')
    with open(mapfile) as fp:
        for line in fp:
            site = line.split()[3]
            chr = line.split()[0]
            filout.write(chr + '\t' + site + '\n')
    filout.close()




for k,v in arrayset.items():
    print(k)
    name = 'phased_' +  k.lstrip("../dataprepper/testset_data/").rstrip('.fam') + '.keepvcf'
    print(name)
    filout = open(name,'w')
    for i in v:
        filout.write(''.join([i,'\n']))
    filout.close()
