import glob
import sys
import os
output=sys.argv[-1]
inputs= sys.argv[1:-1]
if len(inputs) <2:
	raise ValueError('Please give an input')
print(inputs)

curd = os.getcwd() 

names = [curd + '/' + i.rstrip('.bed') for i in inputs]
print(names)

with open(output, 'w+') as f:	#Open allfiles.txt that will be used in merge step
	for t in names:
		f.write(t + '\n')
