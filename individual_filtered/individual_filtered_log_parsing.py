#!/usr/bin/env python 


import sys
from sys import argv
import re

# Takes three arguments at commandline (ensure that paths to each is designated)
script, infile, outfile = argv

# Opens input file from argv
with open (infile, 'r') as logfile:
    log = logfile.read()

action = ['Individual Filtering Step']    
assay = re.findall(r'--bfile ./allele_filtered/([\w]+)', log)#Identifies assay used from input file given to PLINK
assay1 = re.findall(r'--bfile allele_filtered/([\w]+)', log) #Alternative regex for ID of assay
mind = re.findall(r'--(mind\s*[\w.]+)', log)# filter used for individuals with low call rates (individuals missing proportion of variants above this value will be removed)
rem = re.findall(r'([0-9]+) cow removed due to missing genotype data', log)#Extracts number of animals removed (for if single cow is removed)
rem1 = re.findall(r'([0-9]+) cattle removed due to missing genotype data', log)#Like above, but with multiple animals removed
start = re.findall (r'([0-9]+) cattle \(', log)#Identifies starting number of animals loaded for proportion calculation
                
#List of strings from above regexpressions
stats = action + assay + assay1 + mind + rem + rem1 + start 

#Calculates proportion of animals removed
proportion = float(stats[3])/float(stats[4])
stats.append(proportion)

#Searches for and appends to stats list the time when plink filter was run
time = re.findall(r'Start time: ([0-9 a-z A-Z : .]+)',log)
stats.append(time)

#Appends resulting stats list to the output csv file designated by argv 
logfile = open (outfile, 'a')
logfile.write(str(stats) + '\n')
logfile.close()
