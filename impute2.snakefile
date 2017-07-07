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

# IMPGENS =['hol_testset.SNP50.788', 'hol_testset.GGPLD.788', 'hol_testset.F250.197', 'hol_testset.HD.197']
# IMPREFS =['hol_testset.F250.197', 'hol_testset.HD.197']
IMPREFS = ['f250','hd']
IMPGENS = ['snp50','ggpld','f250','hd']


targfiles=[] #troy's runner
for sample in IMPGENS:
	for chr in rangedict.keys():
		chunkcounter=-1
		for chunk in rangedict.get(chr):
			chunkcounter = chunkcounter+1
			file = 'impute2_imputed/'+sample+'.chr'+chr+'.'+str(chunkcounter)+'.phased.impute2'
			targfiles.append(file)
#print(len(targfiles))
targfiles=[] #jesses's runner
for sample in IMPGENS:
	for chr in rangedict.keys():
		chunkcounter=-1
		for chunk in rangedict.get(chr):
			chunkcounter = chunkcounter+1
			file = 'impute2_imputed/'+sample+'.run2.chr'+chr+'.'+str(chunkcounter)+'.phased.impute2'
			targfiles.append(file)

rule impute:
	input:
		targ = expand("impute2_chromosome/run{run}/{sample}.chr{chr}.phased.imputed.gen", sample = IMPGENS, run = 2, chr = list(range(29,30)))

include: "phasing.snakefile"
include: "shapeit.snakefile"

#snakemake -s impute2.snakefile --cores 34  &> imputationrunone_02.txt

#Run 1 (haps, sample) = "vcf_to_haps/{assay}.run{run}.chr{chr}.phased.haps"
#Run 1 (hap, legend, sample) = "vcf_to_hap/{assay}.run{run}.chr{chr}.phased.legend"
#Run 2 (hap, legend, sample) = "impute_input/run{run}/{sample}.chr{chr}.phased.legend"
#Run 2 (haps, sample) = "eagle_phased_assays/run{run}/{sample}.chr{chr}.phased.haps"



# haps_sample_run = {'1': 'vcf_to_haps','12': 'vcf_to_haps', '2': 'eagle_phased_assays','13':'eagle_phased_assays', '4':'shapeit_phased_assays','6':'vcf_to_haps','7':'eagle_phased_assays', '9':'shapeit_phased_assays'}
# haplegendsample_run = {'1':'vcf_to_hap','12': 'vcf_to_hap', '2':'impute_input','13':'impute_input','4':'shapeit_phased_assays/impute_input','6':'vcf_to_hap', '7':'impute_input', '9':'shapeit_phased_assays/impute_input'}

def haps_runlocator(shoein):
	haps_sample_run = {'1': 'vcf_to_haps','12': 'vcf_to_haps', '2': 'eagle_phased_assays','13':'eagle_phased_assays', '4':'shapeit_phased_assays','6':'vcf_to_haps','7':'eagle_phased_assays', '9':'shapeit_phased_assays','14':'shapeit_phased_assays'}
	loc = []
	t = shoein.run
	samp = shoein.sample
	chrom = shoein.chr
	location = haps_sample_run[t] + '/run' + t + '/'+ samp  +'.chr' + chrom + '.phased.haps'
	#print(location)
	return location

def hap_runlocator(shoein):
	haplegendsample_run = {'1':'vcf_to_hap','12': 'vcf_to_hap', '2':'impute_input','13':'impute_input','4':'shapeit_phased_assays/impute_input','6':'vcf_to_hap', '7':'impute_input', '9':'shapeit_phased_assays/impute_input','14':'shapeit_phased_assays/impute_input'}
	loc = []
	for xx in IMPREFS:
		t = shoein.run
		samp = xx
		chrom = shoein.chr
		location = haplegendsample_run[t] +'/run' + t + '/' + samp + '.chr' + chrom + '.phased.haplotypes'
		loc.append(location)
	return loc

def legend_runlocator(shoein):
	haplegendsample_run = {'1':'vcf_to_hap','12': 'vcf_to_hap', '2':'impute_input','13':'impute_input','4':'shapeit_phased_assays/impute_input','6':'vcf_to_hap', '7':'impute_input', '9':'shapeit_phased_assays/impute_input','14':'shapeit_phased_assays/impute_input'}
	loc = []
	for xx in IMPREFS:
		t = shoein.run
		samp = xx
		chrom = shoein.chr
		location = haplegendsample_run[t] +'/run' + t + '/' + samp +'.chr' + chrom + '.phased.legend'
		loc.append(location)
	return loc

rule impute2_refpanel:
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

rule run_impute2_run2: #for parralel phasing
	input:
		hap = "impute2_refpanel/run{run}/merged_refpanel.chr{chr}.{chunk}.phased.hap",	#How are we supposed to expand over a function? Not sure if we can make this as smart as we want it to be?
		legend = "impute2_refpanel/run{run}/merged_refpanel.chr{chr}.{chunk}.phased.legend",
		knownhaps = haps_runlocator,
		maps="impute_maps/imputemap.chr{chr}.map"
	params:
		chunk= chrchunker
	log:
		"logs/impute2/run{run}/{sample}.chr{chr}.{chunk}.log"
	benchmark:
		"benchmarks/impute2/run{run}/{sample}.chr{chr}.{chunk}.benchmark.txt"
	output:
		imputed="impute2_imputed/run{run}/{sample}.chr{chr}.{chunk}.phased.impute2",
		summary="impute2_imputed/run{run}/{sample}.chr{chr}.{chunk}.phased.impute2_summary"
	shell:
		"(impute2 -merge_ref_panels -m {input.maps} -h {input.hap} -l {input.legend} -use_prephased_g -known_haps_g {input.knownhaps} -int {params.chunk} -Ne 200 -o {output.imputed}) > {log}"

# rule run_impute2_round2: #for parralel phasing
# 	input:
# 		hap=expand("vcf_to_hap/{ref}.run{{run}}.chr{{chr}}.hap.gz", ref=IMPREFS),
# 		legend=expand("vcf_to_hap/{ref}.run{{run}}.chr{{chr}}.legend.gz", ref=IMPREFS),
# 		knownhaps="vcf_to_haps/run{run}/{sample}.chr{chr}.hap.gz",
# 		maps="impute_maps/imputemap.chr{chr}.map"
# 	params:
# 		chunk= chrchunker
# 	log:
# 		"logs/impute2/{sample}.{run}.chr{chr}.{chunk}.log"
# 	benchmark:
# 		"benchmarks/impute2/{sample}.{run}.chr{chr}.{chunk}.benchmark.txt"
# 	output:
# 		imputed="impute2_imputed/run{run}/{sample}.chr{chr}.{chunk}.phased.impute2",
# 		summary="impute2_imputed/run{run}/{sample}.chr{chr}.{chunk}.phased.impute2_summary"
# 	shell:
# 		"(impute2 -merge_ref_panels -m {input.maps} -h {input.hap} -l {input.legend} -use_prephased_g -known_haps_g {input.knownhaps} -int {params.chunk} -Ne 200 -o {output.imputed}) > {log}"
# 		#-fill_holes


#print(sample)
#print('sample\n')
rundict = {}
for Run in range(15):
	run= str(Run)
	sampledict = {}
	for sample in IMPGENS:
		filedict = {}
		for chr in rangedict.keys():
			chunkcounter=-1
			flist = []
			for chunk in rangedict.get(chr):
				chunkcounter = chunkcounter+1
				file = 'impute2_imputed/run'+ run + '/'+sample + '.chr'+chr+'.'+str(chunkcounter)+'.phased.impute2' #need to edit this to accept run as a wildcard
				flist.append(file)
				filedict[chr]=flist
		sampledict[sample] = filedict
	rundict[run] = sampledict

def chrfiles(chrom):
	# if chrom.chr == '26':
	# 	print(rundict[chrom.run][chrom.sample][chrom.chr][:4])
	return rundict[chrom.run][chrom.sample][chrom.chr]

#print(rundict['2']['snp50']['26'][:10])

rule cat_chunks:
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
		cat = temp("impute2_chromosome/run{run}/{sample}.chr{chr}.phased.imputed.gen")
	shell:
		"(cat {input.chunks} > {output.cat}) > {log}"#; cp eagle_phased_assays/*.run2*.sample impute2_chromosome/"

#maybe needed?
 #bcftools query -l ./vcf_to_assays/hol_testset.F250.197.1.chr10.vcf > ./vcf_to_assays/hol_testset.F250.197.1.samples

#perl ./bin/vcf2impute_legend_haps.pl -vcf vcf_to_assays/hol_testset.GGPLD.788.1.chr25.vcf -leghap vcf_to_haps/hol_testset.GGPLD.788.1.chr25 -chr 25
