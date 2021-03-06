#!/usr/bin/env python


import sys
from sys import argv
import re

#Takes three arguments on command line to run.  Filepath needed to both input and output files (in different directories).
#Must be a plink .log file in order for regular expressions to work
script, infile, outfile = argv


with open (infile, 'r') as logfile:
    log = logfile.read()

infile= infile.strip('allele_filtered/') #These three seubsequent steps take the input file's name and separate date, #SNPs, #individuals, and Assay to be reported on the top line for calculations, etc.
infile = infile.strip('.log')
metadata = infile.split('.')


action = ['Allele Filtering Step']
#Regular expressions look for metadata in log files.      
geno = re.findall(r'--(geno\s*[\w.]+)', log) #The filter used (loci with proportion of genotypes missing > geno will be removed) 
rem = re.findall(r'([0-9]+) variants removed', log) #This extracts the number of variants removed using the given geno filter
start = re.findall(r'([\w]+) variants, ', log)#Number of variants loaded from input file (for calculating proportion lost)
space = [' ']
time = re.findall(r'Start time: ([0-9 a-z A-Z : .]+)',log)

stats = action + geno + rem + space + start + space #Concatenates lists from previous 4 regular expressions
stats.insert(0, metadata[3])

proportion = int(stats[3])/int(stats[5])#Calculates a proportion of of variants removed by filter.
x = []
x.append(str(proportion))
stats = stats + x + space + time

stats = ', '.join(stats)

logfile = open (outfile, 'a+') #Appends the list generated in previous steps to designated csv outfile
logfile.write(str(stats) + '\n')
logfile.close()
