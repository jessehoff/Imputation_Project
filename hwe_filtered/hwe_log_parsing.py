#!/user/bin/env python


from datetime import datetime
import sys
from sys import argv
import re
script, infile = argv

with open (infile, 'r') as logfile:
    log = logfile.read()

assay = re.findall(r'individual_filtered/([\w_]+)\n', log)
hwe = re.findall(r'--(hwe\s*[\w.]+)', log)
rem = re.findall(r'([0-9]+) variants removed', log)
start = re.findall(r'([\w]+) variants loaded from ', log)
time = re.findall(r'Start time: ([0-9 a-z A-Z : .]+)',log)
stats = assay + hwe + rem + start
percent = int(stats[2])/int(stats[3])
stats.append(percent)
time = re.findall(r'Start time: ([0-9 a-z A-Z : .]+)',log)
stats = stats +time

logfile = open ('hwe_filter_logging.csv', 'a')
logfile.write(str(stats) + '\n')
logfile.close()
