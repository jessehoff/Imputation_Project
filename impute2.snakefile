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

IMPGENS =['hol_testset.SNP50.788', 'hol_testset.GGPLD.788', 'hol_testset.F250.197', 'hol_testset.HD.197']
IMPREFS =['hol_testset.F250.197', 'hol_testset.HD.197']
Phaseruns = {'1': 'vcf_to_haps'}


targfiles=[] #troy's runner
for sample in IMPGENS:
	for chr in rangedict.keys():
		chunkcounter=-1
		for chunk in rangedict.get(chr):
			chunkcounter = chunkcounter+1
			file = 'impute2_imputed/'+sample+'.chr'+chr+'.'+str(chunkcounter)+'.phased.impute2'
			targfiles.append(file)
print(len(targfiles))
targfiles=[] #jesses's runner
for sample in IMPGENS:
	for chr in rangedict.keys():
		chunkcounter=-1
		for chunk in rangedict.get(chr):
			chunkcounter = chunkcounter+1
			file = 'impute2_imputed/'+sample+'.run1.chr'+chr+'.'+str(chunkcounter)+'.phased.impute2'
			targfiles.append(file)
print(len(targfiles))
print(targfiles[:10])

rule targ:
	input:
		#targ=targfiles #This creates the f
		#targ= expand("impute_input/{sample}.chr{chr}.phased.haplotypes", sample = IMPREFS,chr=list(range(20,30)))
		targ = expand("impute2_chromosome/{sample}.run1.chr{chr}.phased.imputed.gen", sample = IMPGENS, chr = list(range(1,30))) #troys seperate phased runner

		#targ=targfiles[:10]

#snakemake -s impute2.snakefile --cores 34  &> imputationrunone_02.txt


def haps_runlocator(shoein):
		t = shoein.run
		samp = shoein.sample
		cra = shoein.chr
		location = Phaseruns[t] + '/' + samp + '.chr' + chr + '.phased.haps'
		return location
def sample_runlocator(shoein):
		t = shoein.run
		samp = shoein.sample
		cra = shoein.chr
		location = Phaseruns[t] + '/' + samp + '.chr' + chr + '.phased.haps'
		return location


rule hap_leg:
	input:
		haps = "eagle_phased_assays/{sample}.chr{chr}.phased.haps",
		sample = "eagle_phased_assays/{sample}.chr{chr}.phased.sample"
	params:
		inprefix = "eagle_phased_assays/{sample}.chr{chr}.phased",
		oprefix = "impute_input/{sample}.chr{chr}.phased"
	log:
		"snake_logs/hap_leg/{sample}.chr{chr}.phased.log"
	benchmark:
		"filter_benchmarks/hap_leg/{sample}.chr{chr}.phased.benchmark.txt"
	output:
		hap = "impute_input/{sample}.{run}.chr{chr}.phased.haplotypes",
		leg = "impute_input/{sample}.{run}.chr{chr}.phased.legend",
		log = "impute_input/logs/{sample}.{run}.chr{chr}.phased.log"
	shell:
		"(shapeit -convert --input-haps {params.inprefix} --output-log {output.log} --output-ref {params.oprefix}) > {log}"




rule run_impute2_round1:
	input:
		haps=expand("impute_input/{sample}.chr{{chr}}.phased.haplotypes", sample=IMPREFS),
		legend=expand("impute_input/{sample}.chr{{chr}}.phased.legend", sample=IMPREFS),
		gens="eagle_phased_assays/{sample}.chr{chr}.phased.haps", #Important to remember that here, we are just referencing the .haps.gz files that are outputted by Eagle directly, instead of putting them through the impute scripts to create .gen.gz files like we might otherwise if not formatted correctly.  This should work
		maps="impute_maps/imputemap.chr{chr}.map"
	params:
		chunk= chrchunker
	log:
		"logs/impute2/{sample}.chr{chr}.{chunk}.log"
	benchmark:
		"benchmarks/impute2/{sample}.chr{chr}.{chunk}.benchmark.txt"
	output:
		imputed="impute2_imputed/{sample}.chr{chr}.{chunk}.phased.impute",
		summary="impute2_imputed/{sample}.chr{chr}.{chunk}.phased.impute2_summary"
	shell:
		"(impute2 -merge_ref_panels -m {input.maps} -h {input.haps} -l {input.legend} -use_prephased_g -known_haps_g {input.gens} -int {params.chunk} -Ne 200 -o {output.imputed}) > {log}"

rule run_impute2_round2: #for parralel phasing
	input:
		hap=expand("vcf_to_hap/{ref}.{{run}}.chr{{chr}}.hap.gz", ref=IMPREFS),
		legend=expand("vcf_to_hap/{ref}.{{run}}.chr{{chr}}.legend.gz", ref=IMPREFS),
		knownhaps="vcf_to_haps/{sample}.{run}.chr{chr}.hap.gz",
		maps="impute_maps/imputemap.chr{chr}.map"
	params:
		chunk= chrchunker
	log:
		"logs/impute2/{sample}.{run}.chr{chr}.{chunk}.log"
	benchmark:
		"benchmarks/impute2/{sample}.{run}.chr{chr}.{chunk}.benchmark.txt"
	output:
		imputed="impute2_imputed/{sample}.run{run}.chr{chr}.{chunk}.phased.impute2",
		summary="impute2_imputed/{sample}.run{run}.chr{chr}.{chunk}.phased.impute2_summary"
	shell:
		"(impute2 -merge_ref_panels -m {input.maps} -h {input.hap} -l {input.legend} -use_prephased_g -known_haps_g {input.knownhaps} -int {params.chunk} -Ne 200 -o {output.imputed}) > {log}"
		#-fill_holes


print(sample)
print('sample\n')
sampledict = {}
for sample in IMPGENS:
	filedict = {}
	for chr in rangedict.keys():
		chunkcounter=-1
		flist = []
		for chunk in rangedict.get(chr):
			chunkcounter = chunkcounter+1
			file = 'impute2_imputed/'+sample+'.run1.chr'+chr+'.'+str(chunkcounter)+'.phased.impute2'
			flist.append(file)
			filedict[chr]=flist
	sampledict[sample] = filedict

def chrfiles(chrom):
	return sampledict[chrom.sample][chrom.chr]

rule cat_chunks:
	input:
		chunks = chrfiles
		#chunks = "impute2_imputed/{sample}.chr{chr}.1.phased.impute2"
	params:
		#star = "impute2_imputed/{sample}.chr{chr}.*.phased.impute2" What is wrong wit this appraoch
	log:
		"logs/cat_chunks/{sample}.run1.chr{chr}.log"
	benchmark:
		"benchmarks/cat_chunks/{sample}.run1.chr{chr}.benchmark.txt"
	output:
		cat = "impute2_chromosome/{sample}.run1.chr{chr}.phased.imputed.gen"
	shell:
		"(cat {input.chunks} > {output.cat}) > {log};"# cp eagle_phased_assays/*.sample impute2_chromosome/"


#maybe needed?
 #bcftools query -l ./vcf_to_assays/hol_testset.F250.197.1.chr10.vcf > ./vcf_to_assays/hol_testset.F250.197.1.samples

#perl ./bin/vcf2impute_legend_haps.pl -vcf vcf_to_assays/hol_testset.GGPLD.788.1.chr25.vcf -leghap vcf_to_haps/hol_testset.GGPLD.788.1.chr25 -chr 25
