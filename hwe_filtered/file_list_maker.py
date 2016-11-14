mport glob
import pandas as pd
import numpy as np

bed = []	#Set up empty lists for each needed file type     
bim = []
fam = []

for i in glob.glob('*.bed'): #Each of these searches for any file in directory that ends in "bed", "bim", or "fam"
    bed.append(i)
for j in glob.glob('*.bim'):
    bim.append(j)
for k in glob.glob('*.fam'):
    fam.append(k)
    
with open('allfiles.txt', 'w+'):	#Open allfiles.txt that will be used in merge step
    data = np.array([bed, bim, fam])	# Reads three lists into a numpy array, then transposes so that assays are rows, and file types are columns
    data = data.T
np.savetxt('allfiles.txt', data, fmt=['%s', '%s','%s']) 	#Writes to text file
