#!/usr/bin/env python  
import sys 
from sys import argv 
import numpy as np
   
script, infile, geno = argv  

with open (infile,'r') as x:     
    frq_input = x.read()  

split = frq_input.split() 
split[0:11] = [] 
nchrobs = split[0::6] 
nchrobs = list(map(int, nchrobs)) 
nchrobs = np.array(nchrobs) 
myint = nchrobs.max() 
callrate = nchrobs/myint 
filtered = 0 

for x in callrate:     
    if x <=1-float(geno):         
        filtered = filtered +1     
    ratio = (filtered/len(callrate))  

if ratio <= 0.05:    
    with open ('test.csv', 'a') as logfile:
        logfile.write(infile + '\t' + geno + '\t' + str(ratio) + '\n')    

else:    
    with open('test.csv', 'a') as logfile:     
        logfile.write(infile + '\t' + geno + '\t' + "This filter removes too many variants!" + '\n')              
    print ("This filter removes too many variants")
