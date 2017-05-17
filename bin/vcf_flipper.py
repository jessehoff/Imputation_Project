import glob as glob
from collections import defaultdict
from collections import Counter
import sys

#Specify chromosome as command line arg
chr = sys.argv[1]

#
allfiles = defaultdict(dict) #initiates dictionary
count = 0
wild = 'assay_vcfs/*.chr' + chr + '.phased.vcf' #identifies the global list of files for the specified chromosome.  Need a more specific identifier here? ex. *chr*.vcf
print(glob.glob(wild))

for file in glob.glob(wild): #iterates through global list of files that are
    fileo = open(file,'r')
    print(file)
    for i in fileo:
        if not i.startswith('#'): #Goes to first line that doesn't begin with a #
            count +=1 #Adds to the counter
            fill =i.split('\t')
            site = fill[1] #stores site number
            order = ','.join(fill[3:5]) #stores Ref/Alt in order
            allfiles[site][file] = order #Dict of dict:  {site:assay:refalt_order}
print(count) #Prints the number of sites that are on multiple assays

diffs = []
for sites,files in allfiles.items(): #iterates through items of allfiles (just sites and assays)
    if len(files) > 1: #if a site is on multiple assays, proceed
        temp = []
        for i in files:
            temp.append(allfiles[sites][i]) #appends Ref/Alt for site from each of the assays that it is a part of
        if len(set(temp))>1: #if there are multiple combinations of ref/alt, append the sites to the diffs list
            diffs.append(sites)
print(len(diffs))	#print the number of ref/alt disagreements

replace = defaultdict(dict)
for i in diffs:
    tally = Counter(allfiles[i].values()).most_common() #For each site that has differeing ref/alt alleles, find the one that is the most common amongst assays, store as tally
    if tally[1][1] == 1: #If the count of the less common refalt ==1
        search = tally[1][0]#Find the less common refalt
        for f,t in allfiles[i].items():
            if t == search:#for every assay for every site in allfiles, if the refalt pair matches the less common refalt found above, then add this entry to the replace dictionary with the more common refalt pair
                replace[f][i]= (tally[0][0],t)#Add this to the dictionary of replacements.  Dict structure: {site:assay:morecommon_ref/alt}
    if tally[1][1] != 1:
        search = tally[1][0]
        for f,t in allfiles[i].items():
            if t == search:
                replace[f][i]= (tally[0][0],t)
flip = {'1|0':'0|1','0|0':'1|1','0|1':'1|0','1|1':'0|0'} #possible combinations of flipping
for i,j in replace.items(): #iterate through replace dict
    fileo = open(i,'r+') #opens outfile
    newname = i.replace('phased','phased.replace')#replace .phased with .phased.replace in new filename
    out = open(newname,'w')
    for p in fileo:
        if not p.startswith('#'): #Go to first non "#" place in file
            fill =p.strip().split('\t')
            site = fill[1]
            if site in j.keys(): #if site is in replace dict
                new = fill[:3] + j[site][0].split(',') + fill[5:9] #replace the refalt order
                tail = [flip[x] for x in fill[9:]]#Flip all of the genotypes in the row from the flip dict above
                jib = new + tail + ['\n'] #add new line after
                out.write('\t'.join(jib)) #write new file
            else:
                out.write(p)
        else:
            out.write(p)
    out.close()
    fileo.close()
