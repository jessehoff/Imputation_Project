#!/usr/bin/env python

import numpy as np
import pandas as pd
from pandas import DataFrame
from pandas import Series
from numpy import nan
from sys import argv
script, infile, outfile = argv

# Takes PLINK ID file as an input.  Splits it into columns, we use only the column with ID in it.
chip = pd.read_table(infile, delim_whitespace=True,usecols=[1], header=None)
chip['is_duplicated'] = chip.duplicated()
duplicates = chip.loc[chip['is_duplicated']==True]
np.savetxt(outfile, duplicates[1].values, fmt='%s')
#with open(outfile, 'w') as file:
