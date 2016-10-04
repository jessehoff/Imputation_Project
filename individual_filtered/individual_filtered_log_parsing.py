#!/usr/bin/env python 

import sys
from sys import argv
script, infile = argv
import re

script, infile = argv

with open (infile, 'r') as logfile:
    log = logfile.read()
    
assay = re.findall(r'--bfile ./allele_filtered/([\w_]+)', log)
mind = re.findall(r'--(mind\s*[\w.]+)', log)
rem = re.findall(r'([0-9]+) cow removed due to missing genotype data', log)
rem1 = re.findall(r'([0-9]+) cattle removed due to missing genotype data', log)
start = re.findall (r'([0-9]+) cattle \(', log)
                

stats = assay + mind + rem + rem1 + start

percent = int(stats[2])/int(stats[3])
stats.append(percent)


logfile = open ('individual_filtered_log.csv', 'a')
logfile.write(str(stats) + '\n')
logfile.close()
