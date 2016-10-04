#!/usr/bin/env python

import sys
from sys import argv
import re
script, infile = argv

with open (infile, 'r') as logfile:
    log = logfile.read()
    
assay = re.findall(r'100.test_([\w_]+)', log)
geno = re.findall(r'--(geno\s*[\w.]+)', log)
rem = re.findall(r'([0-9]+) variants removed', log)
start = re.findall(r'([\w]+) variants, ', log)


stats = assay + geno + rem + start
percent = int(stats[2])/int(stats[3])
stats.append(percent)

logfile = open ('allele_filtered_log.csv', 'a')
logfile.write(str(stats) + '\n')
logfile.close()
