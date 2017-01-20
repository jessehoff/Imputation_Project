import sys
from sys import argv
import re

#Takes three command line arguments (be sure to designate file paths for these as they're likely in different directories)
script, missexed, sexlog, outfile = argv

#Open input file of animals that PLINK determines are missexed 
with open (missexed, 'r') as missexed:
    count = sum(1 for line in missexed)

with open (sexlog, 'r') as logfile:
    log = logfile.read()

sexlog= sexlog.strip('sex_impute/') #These three seubsequent steps take the input file's name and separate date, #SNPs, #individuals, and Assay to be reported on the top line for calculations, etc.
sexlog = sexlog.strip('.log')
metadata = sexlog.split('.')
	
action = ['Missexed Animals Filtered']
impute_sex = re.findall(r'--(impute-sex ycount)\n', log)
rem = re.findall(r'([0-9]+) cattle \(', log)
space = [' ']
time = re.findall(r'Start time: ([0-9 a-z A-Z : .]+)',log)

stats = action + impute_sex + space + [str(count)] + space + rem + space 

stats.insert(0, metadata[3])


#Calculates proportion of animals removed
proportion = int(stats[4])/int(stats[6])
x=[]
x.append(str(proportion))

stats = stats + x + time
stats = ', '.join(stats)


#Appends stats list to the end of specified csv file
csv = open (outfile, 'a')
csv.write(str(stats) + '\n')
csv.close()
