import pandas as pd
import numpy as np
from pandas import DataFrame
import sys
from sys import argv
import re

#Takes script and PLINK .sexcheck file as two arguments on command line
script, infile, outfile = argv

sexcheck = pd.read_table(infile, delim_whitespace=True)	#Reads into pandas datatable
	
sexcheck = sexcheck[sexcheck.PEDSEX !=0]		#Looks for and saves all non-zero values for PEDSEX (animals with sexes assigned on pedigree		  
sexcheck = sexcheck[sexcheck.SNPSEX !=0]
problems = sexcheck[sexcheck.STATUS == 'PROBLEM']	#Looks for and saves rows from previous that have PROBLEM status
my_list =problems['IID'].tolist()
my_list=list(set(my_list))

with open(outfile,'w') as file:		#Writes to designated txt output file, which will act as a wildcard in snakemake (custom name based on the creation of file in that step) -- This is write only, not append, so it will refresh this every time that the step is run
    for item in my_list:
        file.write(item + ' '+item +"\n")
