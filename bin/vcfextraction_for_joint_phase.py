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
    imputeanimal.append(i.split()[0])
samples = glob.glob('../hwe_filtered/*.100*.fam')
samples=sorted(samples, key=lambda sample: int(sample.split('/')[2].split('.')[0]),reverse=True)
samples

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

prin('check *.100*')
maps = glob.glob('../hwe_filtered/*.100*.bim')
maps


for mapfile in maps:
    name = 'phased_' +  mapfile.lstrip("../hwe_filtered/").rstrip('.bim') + '.vcfregion'
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
    name = 'phased_' +  k.lstrip("../hwe_filtered/").rstrip('.fam') + '.keepvcf'
    print(name)
    filout = open(name,'w')
    for i in v:
        filout.write(''.join([i,'\n']))
    filout.close()