#!/usr/bin/env python

from sys import argv
import csv

script, zip_list, samples, ff, outfile = argv

with open (samples, "r") as s:
    samples = s.readlines()
with open (zip_list, "r") as z:
    ziplist = z.readlines()
with open(ff, "r") as f:
    ff = f.readlines()

izip = {}
for i in samples:
    izip[i.split()[0]]=i.split()[1]

pairs = []
for line in ziplist:
    linepair = []
    linepair.append(line.split()[0])
    linepair.append(line.split()[1])
    pairs.append(linepair)

combined = []
for k, v in izip.items():
    for zip, zone in pairs:
        if v in zip:
            combined.append(k)
            combined.append(zone)
combined=[combined[i:i+2] for i in range(0, len(combined), 2)]

fam=[]
for line in ff:
    fam.append(line.split())
for line in fam:
    for ind, code in combined:
        if line[1] == ind:
            line[5]=code


outfile = open (outfile, 'w')
for line in fam:
        entry = ' '.join(line)
        outfile.write(entry+ "\n")
outfile.close()
