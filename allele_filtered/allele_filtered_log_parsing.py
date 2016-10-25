#!/usr/bin/env python


import sys
from sys import argv
import re

#Takes three arguments on command line to run.  Filepath needed to both input and output files (in different directories).
#Must be a plink .log file in order for regular expressions to work
script, infile, outfile = argv


with open (infile, 'r') as logfile:
    log = logfile.read()

action = ['Allele Filtering Step']
#Regular expressions look for metadata in log files.      
assay = re.findall(r'100.test_([\w_]+)', log) #Assay name from input of plink file name
geno = re.findall(r'--(geno\s*[\w.]+)', log) #The filter used (loci with proportion of genotypes missing > geno will be removed) 
rem = re.findall(r'([0-9]+) variants removed', log) #This extracts the number of variants removed using the given geno filter
start = re.findall(r'([\w]+) variants, ', log)#Number of variants loaded from input file (for calculating proportion lost)


stats = action + assay + geno + rem + start #Concatenates strings from previous 4 regular expressions
proportion = int(stats[3])/int(stats[4])#Calculates a proportion of of variants removed by filter.
stats.append(proportion)#Appends to list of strings

time = re.findall(r'Start time: ([0-9 a-z A-Z : .]+)', log) #Regex to grab timestamp from .log file. Appends to list of strings
stats.append(time)


logfile = open (outfile, 'a') #Appends the list generated in previous steps to designated csv outfile
logfile.write( str(stats) + '\n')
logfile.close()
