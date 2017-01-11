#!/usr/bin/env python

from sys import argv
script, id_file, duplicate_file, old_ped, new_ped  = argv
#User inputs are 1) the script   2) ID file with possible duplicates   3) File to write duplicate ID's to   4) ped file with possible duplicates   5) New ped file to write non-duplicated results to.

ids=[]
chip = open(id_file, 'r').readlines()
outfile = open(duplicate_file, 'w')
for line in chip:		
    id = line.split()[1]		#Parses ID file, looking only at column with ID's
    if id not in ids:
        ids.append(id)			#Upon reading an ID, adds to ids list, which means it cannont be written to the outfile
    else: 
        outfile.write(line.split()[1] + '\n')
outfile.close()


dup={}					#Duplicates are added to this dictionary, and assigned the value 0
file = open(duplicate_file, 'r')
for line in file:
    line = line.strip()
    dup.update({line:0})

    
outfile =  open(new_ped, 'w')    
ped = open(old_ped, 'r').readlines()
for x in ped:			
    id = x.split()[0]			# using the IID from ped file, it can be written to the new ped file if it's value in dict is 0
    if id in dup:			
            if dup[id] == 0: 
                outfile.write(x) 
                dup[id] +=1		#upon being written, an ID's value in dict increases to 1, preventing it from being written to the new output
    else:
        outfile.write(x)
outfile.close()
