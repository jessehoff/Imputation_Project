from itertools import tee
def parwise(iterable):
	a,b = tee(iterable)
	next(b, None)
	return zip(a,b)

def chunks(end):
	return [str(i[0]+1)+ ' ' + str(i[1]) for i in parwise(range(0,end,5000000))]

rangedict = {'1':chunks(158322647),'2':chunks(136914030),'3':chunks(121412973), '4':chunks(120786530), '5':chunks(121186724), '6':chunks(119454666), '7':chunks(112628884), '8':chunks(113380773), '9':chunks(105701306), '10':chunks(104301732), '11':chunks(107282960), '12':chunks(91155612), '13':chunks(84230359), '14':chunks(84628243), '15':chunks(85272311), '16':chunks(81720984), '17':chunks(75149392), '18':chunks(65999195), '19':chunks(64044783), '20':chunks(71992748), '21':chunks(71594139), '22':chunks(61379134), '23':chunks(52467978), '24':chunks(62685898), '25':chunks(43879707), '26':chunks(51680135), '27':chunks(45402893), '28':chunks(46267578), '29':chunks(51504286)} #, '30':chunks(143032828)}
#describes the length of each chromsome so taht it can be chunked into equal sized chunks.

def chrchunker(WC):
	return rangedict[WC.chr][int(WC.chunk)]

#include: 'phasing.snakefile'
# IMPGENS =['hol_testset.SNP50.788', 'hol_testset.GGPLD.788', 'hol_testset.F250.197', 'hol_testset.HD.197']
# IMPREFS =['hol_testset.F250.197', 'hol_testset.HD.197']
IMPREFS = ['f250','hd']
IMPGENS = ['snp50','ggpld','f250','hd']

#include: "phasing.snakefile"
#include: "impute2.snakefile"

rule imp4:
	input:
		targ = expand("impute4_chromosome/run31/{sample}.chr{chr}.imputed.gen.gz", sample = IMPGENS, chr = list(range(1,30)))

# def haps_runlocator(WC):
# 	shapeit = ['17', '18', '19']
# 	eagle_assay = ['30']
# 	eagle_combined = []
# 	loc = []
# 	for xx in IMPREFS:
# 		if WC.run in eagle_combined:
# 			directory = 'vcf_to_haps'
# 		if WC.run in eagle_assay:
# 			directory = 'eagle_phased_assays'
# 		if WC.run in shapeit:
# 			directory = 'shapeit_phased_assays'
# 	location = directory + '/run' + WC.run + '/'+ WC.sample  +'.chr' + WC.chr + '.phased.haps'
# 	#print(location)
# 	return location
#
# def hap_runlocator(WC):
# 	shapeit = ['17', '18', '19']
# 	eagle_assay = ['30']
# 	eagle_combined = []
# 	loc = []
# 	for xx in IMPREFS:
# 		if WC.run in eagle_combined:
# 			directory = 'vcf_to_hap'
# 		if WC.run in eagle_assay:
# 			directory = 'impute_input'
# 		if WC.run in shapeit:
# 			directory = 'shapeit_phased_assays/impute_input'
# 		location = directory + '/run' + WC.run  + '/'+ xx  +'.chr' + WC.chr + '.phased.haplotypes'
# 		loc.append(location)
# 	return loc
#
# def legend_runlocator(WC):
# 	shapeit = ['17', '18', '19']
# 	eagle_assay = ['30']
# 	eagle_combined = []
# 	loc = []
# 	for xx in IMPREFS:
# 		if WC.run in eagle_combined:
# 			directory = 'vcf_to_hap'
# 		if WC.run in eagle_assay:
# 			directory = 'impute_input'
# 		if WC.run in shapeit:
# 			directory = 'shapeit_phased_assays/impute_input'
# 		location = directory + '/run' + WC.run  + '/'+ xx  +'.chr' + WC.chr + '.phased.legend'
# 		loc.append(location)
# 	return loc
#
def haps_runlocator(WC):
	haps_sample_run = {'30':{'vcf_to_haps':'1'}, '31':{'eagle_phased_assays':'2'}}
	directory = list(haps_sample_run[WC.run])[0]
	location = directory + '/run' + haps_sample_run[WC.run][directory] + '/'+ WC.sample  +'.chr' + WC.chr + '.phased.haps.gz'
	return location

def hap_runlocator(shoein):
	haplegendsample_run = {'30':'vcf_to_hap/run1/','31': 'impute_input/run2/'}
	loc = []
	for xx in IMPREFS:
		t = shoein.run
		samp = xx
		chrom = shoein.chr
		location = haplegendsample_run[t] + samp + '.chr' + chrom + '.phased.haplotypes'
		loc.append(location)
	return loc

def legend_runlocator(shoein):
	haplegendsample_run = {'30':'vcf_to_hap/run1/','31': 'impute_input/run2/'}
	loc = []
	for xx in IMPREFS:
		t = shoein.run
		samp = xx
		chrom = shoein.chr
		location = haplegendsample_run[t] + samp +'.chr' + chrom + '.phased.legend'
		loc.append(location)
	return loc


rule impute4_refpanel:
	input:
		hap = hap_runlocator,
		legend = legend_runlocator,
		maps="impute_maps/imputemap.chr{chr}.map"
	log:
		"logs/refpanel_impute/run{run}/merged_refpanel.chr{chr}.phased.log"
	params:
		chunk = chrchunker,
		oprefix = "impute2_refpanel/run{run}/merged_refpanel.chr{chr}.{chunk}.phased"
	benchmark:
		"benchmarks/impute2_refpanel/run{run}/merged_refpanel.chr{chr}.{chunk}.phased.benchmark.txt"
	output:
		refhap = "impute2_refpanel/run{run}/merged_refpanel.chr{chr}.{chunk}.phased.hap",
		refleg = "impute2_refpanel/run{run}/merged_refpanel.chr{chr}.{chunk}.phased.legend"
	shell:
		"(impute2 -merge_ref_panels_output_ref {params.oprefix} -m {input.maps} -h {input.hap} -l {input.legend} -int {params.chunk} -Ne 200 -o {params.oprefix}) > {log}"

rule run_impute4: #for parralel phasing
	input:
		hap = "impute2_refpanel/run{run}/merged_refpanel.chr{chr}.{chunk}.phased.hap",	#How are we supposed to expand over a function? Not sure if we can make this as smart as we want it to be?
		legend = "impute2_refpanel/run{run}/merged_refpanel.chr{chr}.{chunk}.phased.legend",
		knownhaps = haps_runlocator,
		maps="impute_maps/imputemap.chr{chr}.map"
	params:
		chunk= chrchunker,
		oprefix="impute4_imputed/run{run}/{sample}.chr{chr}.{chunk}"
	log:
		"logs/impute4/run{run}/{sample}.chr{chr}.{chunk}.log"
	benchmark:
		"benchmarks/impute4/run{run}/{sample}.chr{chr}.{chunk}.benchmark.txt"
	output:
		imputed="impute4_imputed/run{run}/{sample}.chr{chr}.{chunk}.gen.gz",
		#summary="impute2_imputed/run{run}/{sample}.chr{chr}.{chunk}.phased.impute2_summary"
	shell:
		"(impute4.r265.1 -m {input.maps} -h {input.hap} -l {input.legend} -g {input.knownhaps} -int {params.chunk} -Ne 200 -no_maf_align -o_gz -o {params.oprefix}) > {log}"

rundict = {}
for Run in range(50):
	run= str(Run)
	sampledict = {}
	for sample in IMPGENS:
		filedict = {}
		for chr in rangedict.keys():
			chunkcounter=-1
			flist = []
			for chunk in rangedict.get(chr):
				chunkcounter = chunkcounter+1
				file = 'impute4_imputed/run'+ run + '/'+sample + '.chr'+chr+'.'+str(chunkcounter)+'.gen.gz' #need to edit this to accept run as a wildcard
				flist.append(file)
				filedict[chr]=flist
		sampledict[sample] = filedict
	rundict[run] = sampledict

def chrfiles(chrom):
	# if chrom.chr == '26':
	# 	print(rundict[chrom.run][chrom.sample][chrom.chr][:4])
	return rundict[chrom.run][chrom.sample][chrom.chr]

rule impute4_cat_chunks:
	input:
		chunks = chrfiles
		#chunks = "impute2_imputed/{sample}.chr{chr}.1.phased.impute2"
	params:
		#star = "impute2_imputed/{sample}.chr{chr}.*.phased.impute2" What is wrong wit this appraoch
	log:
		"logs/cat_chunks/run{run}/{sample}.chr{chr}.log"
	benchmark:
		"benchmarks/cat_chunks/run{run}/{sample}.chr{chr}.benchmark.txt"
	output:
		cat = temp("impute4_chromosome/run{run}/{sample}.chr{chr}.imputed.gen.gz")
	shell:
		"(cat {input.chunks} > {output.cat}) > {log}"
