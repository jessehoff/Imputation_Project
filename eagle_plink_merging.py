import glob
import pandas as pd
import numpy as np
import sys
import os
output=sys.argv[1]

curd = os.getcwd() 

bed = []	#Set up empty lists for each needed file type     
bim = []
fam = []

for i in glob.glob(curd  + '/eagle_plink/*.bed'): #Each of these searches for any file in directory that ends in "bed", "bim", or "fam"
    bed.append(i.rstrip('.phased'))

with open(output, 'w+') as f:	#Open allfiles.txt that will be used in merge step
   for t in bed:
       f.write(t + '\n')
