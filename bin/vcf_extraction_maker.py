from sys import argv

script, bim, fam, snps, inds = argv

with open(bim, 'r') as b:
    stringlist = [line.split()[0] + '\t' + line.split()[3] for line in b]
with open(snps, 'w') as outfile:
    [outfile.write(x + '\n') for x in stringlist]

with open(fam, 'r') as f:
    famlist = [line.split()[1] for line in f]
with open(inds, 'w') as outfile:
    [outfile.write(x + '\n') for x in famlist]
