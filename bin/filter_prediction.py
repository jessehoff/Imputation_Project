#!/usr/bin/env python

import sys
from sys import argv
import numpy as np

script, infile, geno = argv


with open (infile,'r') as x:
    frq_input = x.read()
    
split = frq_input.split()
split[0:11] = []
nchrobs = split[0::6]
nchrobs = list(map(int, nchrobs))
nchrobs = np.array(nchrobs)
myint=200
callrate = nchrobs/myint
callrate
filtered = 0
for x in callrate:
    if x <=1-float(geno):
        filtered = filtered +1
ratio = (filtered/len(callrate))

if ratio <= 0.05:
  print (ratio)
else: 
  print ("This filter removes too many variants")
