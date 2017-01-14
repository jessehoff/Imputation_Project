#!/user/bin/env python

import sys
from sys import argv
import re

#Takes three command line arguments (be sure to designate file paths for these as they're likely in different directories)
script, infile, outfile = argv

#Open input log file from PLINK filter
with open (infile, 'r') as logfile:
    log = logfile.read()

action = ['HWE Filtering Step']
#Regular expressions for locating metadata in log file
#assay = re.findall(r'individual_filtered/([\w_]+)\n', log) #Locates assay ID
hwe = re.findall(r'--(hwe\s*[\w.]+)\n', log)#Value of hwe filter
rem = re.findall(r'([0-9]+) variants removed', log)#Number of variants removed by filter
start = re.findall(r'([\w]+) variants loaded from ', log)#Number of variants loaded before filtering

#Concatenates regular expressions above into "stats" list 
stats = action + hwe + rem + start

#Calculates proportion of variants removed compared to total then appends to "stats"
percent = int(stats[2])/int(stats[3])
stats.append(percent)

#Locates time when filter was run by PLINK then appends this to stats list
time = re.findall(r'Start time: ([0-9 a-z A-Z : .]+)',log)
stats = stats +time

#Appends stats list to the end of specified csv file
logfile = open (outfile, 'a')
logfile.write(str(stats) + '\n')
logfile.close()
