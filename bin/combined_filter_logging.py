#!/usr/bin/env python

import re as re
import csv
from sys import argv

script, duplicates, snps, individuals, hwe, sex, missexed, outfile = argv

with open(duplicates, 'r') as d:
	dup = d.readlines()
data = duplicates.split('/')[-1].strip('.txt')
start = data.split('.')
beforelist = '\t'.join([data, 'Before Filtering Steps','','',start[0],'','', start[2],'', start[1]])
duplist = '\t'.join([data, 'Duplicate Individuals Removed','','','','',str(len(dup)), start[2], str(len(dup)/int(start[2]))])

with open(snps, 'r') as s:
	snp = s.read()
geno = re.findall(r'--(geno\s*[\w.]+)', snp)[0]
rem = re.findall(r'([0-9]+) variants removed', snp)[0]
startsnps = re.findall(r'([\w]+) variants, ', snp)[0]
time = re.findall(r'Start time: ([0-9 a-z A-Z : .]+)',snp)[0]
snplist = '\t'.join([data, 'SNP Call Rate', geno, rem, startsnps, str(int(rem)/int(startsnps)),'','','', time])

with open(individuals, 'r') as i:
	ind = i.read()
mind = re.findall(r'--(mind\s*[\w.]+)', ind)[0]
rem = re.findall(r'([0-9]+)\s+[\w.]+ removed due to missing genotype data', ind)[0]
start = re.findall (r'([0-9]+) cattle \(', ind)[0]
time = re.findall(r'Start time: ([0-9 a-z A-Z : .]+)',ind)[0]
indlist = '\t'.join([data, 'Individual Call Rate',mind ,'','','',rem , start, str(int(rem)/int(start)), time])

with open(hwe, 'r') as h:
	hardy = h.read()
hw = re.findall(r'--(hwe\s+\S+)\n', hardy)[0]
rem = re.findall(r'([0-9]+) variants removed', hardy)[0]
start = re.findall(r'([\w]+) variants loaded from ', hardy)[0]
time = re.findall(r'Start time: ([0-9 a-z A-Z : .]+)',hardy)[0]
hwelist = '\t'.join([data, 'Hardy-Weinberg', hw, rem, start, str(int(rem)/int(startsnps)),'','','', time])

with open(sex, 'r') as s: #This is the sexcheck file that PLINK generates
	ms = s.read()
with open(missexed, 'r') as m: #This is the file that has the IDs of missexed animals that my script creates
	rem = str(len(m.readlines()))
impute_sex = re.findall(r'--(impute-sex ycount)\n', ms)[0]
start = re.findall(r'([0-9]+) cattle \(', ms)[0]
time = re.findall(r'Start time: ([0-9 a-z A-Z : .]+)',ms)[0]
missexedlist = '\t'.join([data, 'Missexed Individuals',impute_sex ,'','','',rem , start, str(int(rem)/int(start)), time])

with open(outfile, 'a') as test:
    test.write(beforelist+'\n'+duplist+'\n'+snplist+'\n'+indlist+'\n'+hwelist+'\n'+missexedlist+'\n')
