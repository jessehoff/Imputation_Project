#!/usr/bin/env python
import sys
from sys import argv
import numpy as np

script, infile, mind = argv

with open(infile, 'r') as x:
    imiss = x.read()
split =imiss.split()
split[0:11] = []
fmiss = split[0::6]
fmiss
missing=np.array(fmiss)
filtered = 0

for x in missing:
    if float(x) >= float(mind):
        filtered = filtered + 1
    
if filtered/len(missing) > 0.01:
    with open('predictions.csv', 'a') as log:
        log.write(infile + '\t' + str(mind) + '\t' + "This filter removes too many animals -- " + str(filtered) + '\n')    
    print (infile + '\t' + str(mind) + '\t' + "This filter removes too many animals -- " + str(filtered))
else:
    with open('predictions.csv', 'a') as log:
        log.write(infile + '\t' + str(mind) + '\t' + str(filtered) + " animal(s) removed" + '\n')

    #print (infile + '\t' + str(mind) + '\t' + str(filtered) + " animal(s) removed")
