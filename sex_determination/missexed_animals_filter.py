import pandas as pd
import numpy as np
from pandas import DataFrame
import sys
from sys import argv
import re

#Takes script and PLINK .sexcheck file as two arguments on command line
script, infile = argv

sexcheck = pd.read_table(infile, delim_whitespace=True)	#Reads into pandas datatable
	
sexcheck = sexcheck[sexcheck.PEDSEX !=0]		#Looks for and saves all non-zero values for PEDSEX (animals with sexes assigned on pedigree		  
problems = sexcheck[sexcheck.STATUS == 'PROBLEM']	#Looks for and saves rows from previous that have PROBLEM status
my_list =problems['IID'].tolist()

with open('missexed_animals.txt','w') as file:		#Writes to missexed_animals.txt -- This is write only, not append, so it will refresh this every time that the step is run
    for item in my_list:
        file.write(item + ' '+item +"\n")
