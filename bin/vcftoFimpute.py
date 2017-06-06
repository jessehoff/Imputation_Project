#! /usr/bin/env python
import subprocess
import os
import sys
import numpy as np
import gzip

vcf_dir_file=sys.argv[1]
newend= sys.argv[2]
newstart = sys.argv[4]
colnum = sys.argv[3]
cut = vcf_dir_file.rstrip('vcf') + newend
print(cut) 
put = newstart + cut.split('/')[1]
print(put)
vcf= open(vcf_dir_file,'rb')
out = open(put,'w')




print 'yes'


genskey = {'./.':5, '0/0':'0', '0/1':'1' ,'1/1':'2'}
genslist = []
for i in vcf:
	if i[:6] == '#CHROM':
		ids = i.split()[9:]
		print i[:40]
	if i[0]!="#":	
		gens = [genskey[j] for j in i.strip().split()[9:]]
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



