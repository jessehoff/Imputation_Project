import sys
from sys import argv
import re

#Takes three command line arguments (be sure to designate file paths for these as they're likely in different directories)
script, infile, outfile = argv

#Open input file of animals that PLINK determines are missexed
with open (infile, 'r') as logfile:
    count = sum(1 for line in logfile) #counts lines and reports this as the number of duplicate animals in the dataset

#Steps here hard code in action and file name that will be used for the beginning of the 
action = ['Duplicate Individuals Removed']
starting =['Before Filtering Steps']
space = [' ']
infile= infile.strip('duplicates_filtered/dup_ids/') #These three seubsequent steps take the input file's name and separate date, #SNPs, #individuals, and Assay to be reported on the top line for calculations, etc. 
infile = infile.strip('.txt')
metadata = infile.split('.')

#Establishes list to be written as first line of each file's logging information.Beginning statistics.

start_info =  starting + space + space 
start_info.insert(0, metadata[3])
start_info.append(metadata[0])
start_info = start_info + space + space
start_info.append(metadata[2])
start_info = start_info + space
start_info.append(metadata[1])
start_info = ', '.join(start_info)


removed = []			#Puts the count from above into a "list" so that it can be the same type as other data
removed.append(str(count))

stats =  action + space + space + space + space + removed #Adds the action from above to stats, and reports number of removed entries
stats.append(metadata[2])
stats.insert(0,metadata[3])

proportion = int(stats[6])/int(stats[7])
x=[]
x.append(str(proportion))
stats = stats + x
stats = ', '.join(stats)


#Appends stats list to the end of specified csv file
logfile = open (outfile, 'a')
logfile.write(start_info + '\n' + str(stats) + '\n')
logfile.close()

