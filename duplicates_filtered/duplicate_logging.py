import sys
from sys import argv
import re

#Takes three command line arguments (be sure to designate file paths for these as they're likely in different directories)
script, infile, outfile = argv

#Open input file of animals that PLINK determines are missexed
with open (infile, 'r') as logfile:
    count = sum(1 for line in logfile)

action = ['Duplicate Individuals Removed']

action.append(count)

stats = action
#Appends stats list to the end of specified csv file
logfile = open (outfile, 'a')
logfile.write(str(stats) + '\n')
logfile.close()

