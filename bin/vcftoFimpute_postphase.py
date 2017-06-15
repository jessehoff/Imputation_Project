#! /usr/bin/env python
import subprocess
import os
import sys
import numpy as np
import gzip

vcf_dir_file=sys.argv[1]
newend= sys.argv[2]
colnum = sys.argv[3]
put = vcf_dir_file.rstrip('vcf') + newend 

vcf= open(vcf_dir_file,'rb')
out = open(put,'w')




print 'yes'


genslist = {'./.':5, '0/0':'0', '0/1':'1' ,'1/1':'2'}
for i in vcf:
	if i[:6] == '#CHROM':
		ids = i.split()[9:]
		print i[:40]
	if i[0]!="#":	
		gens = [str(int(i[0])+int(i[2])) for i in i.strip().split()[9:]]
		genslist.append(gens)
print 'lines read'
table = np.array(genslist)
del genslist
table = table.T
print 'done'

table = table.tolist()

for i, j in enumerate(table):
	out.write(''.join([ids[i],'  ' ,colnum ,'  ',''.join([str(p) for p in j]),'\n']))

out.close()



