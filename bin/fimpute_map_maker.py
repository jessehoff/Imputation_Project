#! /usr/bin/env python
import subprocess
import os
import sys

from collections import defaultdict
import gzip

filelist =  sys.argv[1:-2] # a list of bim files for each assay included in the fimpute genotype file
chr = sys.argv[-2]
print chr


sitedict = defaultdict(lambda: ['0','0','0','0']) 
for p,filen in enumerate(filelist):
	hand = open(filen,'r')
	for i, j in enumerate(hand,start=1):
		pos = int(j.split()[3])
		sitedict[pos][p] = str(i)
sites = sitedict.keys()
sites.sort()
print(sites[:10])
print(sitedict[sites[0]])

outn = sys.argv[-1]
out = open(outn,'w')
out.write('snpid, chr ,pos, ch1, ch2, etc \n')


for i in sites:
	nam = str(i)
	name = ':'.join([chr,nam])
	ranks = '\t'.join(sitedict[i])
	out.write('\t'.join([name,chr,nam,ranks,'\n']))

out.close()
