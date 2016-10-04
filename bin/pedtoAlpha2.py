#The genotype file should contain one line for each individual, a first column with the individual’s identifier and as many columns as SNP in the chromosome.
import sys

def twosums(liner):
    gens = []
    for x,y in twogrouped(liner):
        gens.append(str(int(x)+int(y))) #takes two alleles at a time (from an indinvidual) and adds them
    return ','.join(gens)

def twogrouped(iterable):
    return zip(*[iter(iterable)]*2)

fname = sys.argv[1]
outname = sys.argv[2]
file = open(fname)
out = open(outname,'w')
names = []
for line in file:
    total = line.split()
    header = total[:6]
    genos = twosums(total[6:])
    names.append(header[0])
    fin=header[0] + ','+  genos + '\n'
    out.write(fin)
gencount=str(len(genos.split(',')))
print(gencount)
out.close()

#writes a pedfile with null parents ftm
pedf = 'pedtest.txt'
pedfile = open(pedf,'w') # need to not redo every time i run
for an in names:
    pedfile.write(','.join([an,'0,0\n']))
pedfile.close()

specf = open('/home/jlh4df/AlphaImpute1.6_Linux/AlphaImputeSpec.txt')
spec = specf.readlines()
specdict = {}
for t in spec:
    if len(t.split(','))==2:
        k,v = t.split(',')
        specdict[k] = v

specdict['NumberSnp\t\t\t\t'] = gencount + '\n'
specdict['GenotypeFile\t\t\t\t'] = outname + '\n'
specdict['PedigreeFile\t\t\t\t'] = pedf + '\n'

specoutn= outname.replace('.Alpha','.spec')
specout = open(specoutn,'w')

for line in spec:
    name = line.split(',')[0]
    if name in specdict.keys():
        #print(','.join([name,specdict[name]]))
        specout.write(','.join([name,specdict[name]]))
    else:
        specout.write(line)
specout.close()