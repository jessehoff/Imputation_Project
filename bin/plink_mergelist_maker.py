#!/usr/bin/env python

import glob
from sys import argv

script, outlist = argv

ids = glob.glob("plink_imputed/*.chr*imputed.fam")
print(ids)

with open(outlist, 'w') as list:
	for xx in ids:
		prefix = xx.strip(".fam")
		list.write(prefix+'\n')
