#!/usr/bin/env python

from sys import argv
script, id_file, duplicate_file, old_ped, new_ped  = argv

ids=[]
chip = open(id_file, 'r').readlines()
outfile = open(duplicate_file, 'w')
for line in chip:
    id = line.split()[1]
    if id not in ids:
        ids.append(id)
    else: 
        outfile.write(line.split()[1] + '\n')
outfile.close()


dup={}
file = open(duplicate_file, 'r')
for line in file:
    line = line.strip()
    dup.update({line:0})

    
outfile =  open(new_ped, 'w')    
ped = open(old_ped, 'r').readlines()
for x in ped:
    id = x.split()[0]
    if id in dup:
            if dup[id] == 0: 
                outfile.write(x) 
                dup[id] +=1
    else:
        outfile.write(x)
outfile.close()
