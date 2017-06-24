IMPREFS =['hd', 'f250']
SAMPLES = ['hd', 'f250', 'ggpld', 'snp50']

from itertools import tee
def parwise(iterable):
	a,b = tee(iterable)
	next(b, None)
	return zip(a,b)
def chunks(end):
	return [str(i[0]+1)+ ' ' + str(i[1]) for i in parwise(range(0,end,5000000))]
rangedict = {'1':chunks(158322647),'2':chunks(136914030),'3':chunks(121412973), '4':chunks(120786530), '5':chunks(121186724), '6':chunks(119454666), '7':chunks(112628884), '8':chunks(113380773), '9':chunks(105701306), '10':chunks(104301732), '11':chunks(107282960), '12':chunks(91155612), '13':chunks(84230359), '14':chunks(84628243), '15':chunks(85272311), '16':chunks(81720984), '17':chunks(75149392), '18':chunks(65999195), '19':chunks(64044783), '20':chunks(71992748), '21':chunks(71594139), '22':chunks(61379134), '23':chunks(52467978), '24':chunks(62685898), '25':chunks(43879707), '26':chunks(51680135), '27':chunks(45402893), '28':chunks(46267578), '29':chunks(51504286)} #, '30':chunks(143032828)}
def chrchunker(WC):
	return rangedict[WC.chr][int(WC.chunk)]

include: 'phasing.snakefile'
include: 'impute2.snakefile'
rule mm_target:
	input:
		targ = expand("minimac_imputed/run{run}/{sample}.chr{chr}.imputed.dose.vcf.gz", sample = SAMPLES, run = 5, chr = list(range(29,30)))

#hapleg_location = {'impute_input':'2'}

def haps_runlocator(WC):
	haps_sample_run = {'5':{'eagle_phased_assays':'2'}, '10':{'eagle_phased_assays':'7'}, '8':{'vcf_to_haps':'1'}, '11':{'vcf_to_haps':'6'}}
	t = WC.run
	samp = WC.sample
	chrom = WC.chr
	directory = list(haps_sample_run[t])[0]
	location = directory + '/run' + haps_sample_run[t][directory] + '/'+ samp  +'.chr' + chrom + '.phased.haps'
	return location

def haps_runprefix(WC):
	haps_sample_run = {'5':{'eagle_phased_assays':'2'}, '10':{'eagle_phased_assays':'7'}, '8':{'vcf_to_haps':'1'}, '11':{'vcf_to_haps':'6'}}
	t = WC.run
	samp = WC.sample
	chrom = WC.chr
	directory = list(haps_sample_run[t])[0]
	location = directory + '/run' + haps_sample_run[t][directory] + '/'+ samp  +'.chr' + chrom + '.phased'
	return location

def sample_runlocator(WC):
	haps_sample_run = {'5':{'eagle_phased_assays':'2'}, '10':{'eagle_phased_assays':'7'}, '8':{'vcf_to_haps':'1'}, '11':{'vcf_to_haps':'6'}}
	loc = []
	for xx in IMPREFS:
		t = WC.run
		samp = xx
		chrom = WC.chr
		directory = list(haps_sample_run[t])[0]
		location = directory + '/run' + haps_sample_run[t][directory] + '/'+ samp  +'.chr' + chrom + '.phased.sample'
		loc.append(location)
	return loc

def hap_runlocator(WC):
	haplegendsample_run = {'5':'impute2_refpanel/run2/', '8':'impute2_refpanel/run1/', '10':'impute2_refpanel/run7/', '11':'impute2_refpanel/run6/'}
	t = WC.run
	chrom = WC.chr
	chunk = WC.chunk
	location = haplegendsample_run[t] + 'merged_refpanel.chr' + chrom + '.' + chunk + '.phased.hap'
	return location

def legend_runlocator(WC):
	haplegendsample_run = {'5':'impute2_refpanel/run2/', '8':'impute2_refpanel/run1/', '10':'impute2_refpanel/run7/', '11':'impute2_refpanel/run6/'}
	t = WC.run
	chrom = WC.chr
	chunk = WC.chunk
	location = haplegendsample_run[t] + 'merged_refpanel.chr' + chrom + '.' + chunk + '.phased.legend'
	return location

rule reformat_leg:
	input:
		refleg = legend_runlocator
	log:
		"logs/reformat_leg/run{run}/merged_refpanel.chr{chr}.{chunk}.phased.log"
	params:
		chr = "{chr}"
	benchmark:
		"benchmarks/reformat_leg/run{run}/merged_refpanel.chr{chr}.{chunk}.phased.benchmark.txt"
	output:
		newleg = "minimac_ref_panel/ref_panels/run{run}/merged_refpanel.chr{chr}.{chunk}.phased.reformatted.legend"
	shell:
		"(python bin/bcftools_legend_converter.py {params.chr} {input.refleg} {output.newleg}) > {log}"

rule cat_sample:
	input:
		#samples = expand("eagle_phased_assays/{sample}.chr{{chr}}.run{run}.phased.sample", sample = IMPREFS)
		samples = sample_runlocator
	log:
		"logs/reformat_samp/run{run}/merged_refpanel.chr{chr}.phased.log"
	benchmark:
		"benchmarks/reformat_samp/run{run}/merged_refpanel.chr{chr}.phased.benchmark.txt"
	output:
		newsample = "minimac_ref_panel/ref_panels/run{run}/merged_refpanel.chr{chr}.phased.reformatted.sample"
	shell:
		"(python bin/bcftools_sample_converter.py {input.samples} {output.newsample}) > {log}"

rule ref_panel_vcf:
	input:
		hap = hap_runlocator,
		leg = "minimac_ref_panel/ref_panels/run{run}/merged_refpanel.chr{chr}.{chunk}.phased.reformatted.legend",
		samp = "minimac_ref_panel/ref_panels/run{run}/merged_refpanel.chr{chr}.phased.reformatted.sample"
	log:
		"logs/ref_panel_vcf/run{run}/merged_refpanel.chr{chr}.phased.log"
	benchmark:
		"benchmarks/ref_panel_vcf/run{run}/merged_refpanel.chr{chr}.{chunk}.phased.benchmark.txt"
	output:
		vcf = temp("minimac_ref_panel/ref_panel_chunks_vcf/run{run}/merged_refpanel.chr{chr}.{chunk}.phased.vcf")
	shell:
		"(bcftools convert -H {input.hap},{input.leg},{input.samp} -o {output.vcf})"


def chrfiles(chrom):
	rundict = {}
	for Run in range(15):
		run= str(Run)
		filedict = {}
		for chr in rangedict.keys():
			chunkcounter=-1
			flist = []
			for chunk in rangedict.get(chr):
				chunkcounter = chunkcounter+1
				file = 'minimac_ref_panel/ref_panel_chunks_vcf/run'+ run + '/merged_refpanel.chr'+chr+'.'+str(chunkcounter)+'.phased.vcf' #need to edit this to accept run as a wildcard
				flist.append(file)
				filedict[chr]=flist
		rundict[run] = filedict
	return rundict[chrom.run][chrom.chr]

rule cat_vcfs:
	input:
		vcf = chrfiles
	log:
		"logs/cat_vcfs/run{run}/merged_refpanel.chr{chr}.phased.log"
	benchmark:
		"benchmarks/cat_vcfs/run{run}/merged_refpanel.chr{chr}.phased.benchmark.txt"
	output:
		vcf = temp("minimac_ref_panel/ref_panel_chrom_vcf/run{run}/merged_refpanel.chr{chr}.phased.vcf")
	shell:
		"(bcftools concat {input.vcf} -o {output.vcf})"


rule ref_panel_impacc:
	input:
		true = "ref_vcfs/F250_HD_merged.chr{chr}.pickle",
		imputed = "minimac_ref_panel/ref_panel_chrom_vcf/run{run}/merged_refpanel.chr{chr}.phased.vcf",
		frq = "ref_vcfs/F250_HD_merged.chr{chr}.frq",
	params:
		chrom = "{chr}",
		acc = "minimac_ref_panel/ref_panel_impacc/merged_refpanel.txt"
	log:
		"logs/ref_imp_acc/run{run}/merged_refpanel.chr{chr}.txt"
	benchmark:
		#"benchmarks/imp_acc/{sample}.chr{chr}.benchmark.txt"
		"benchmarks/ref_imp_acc/run{run}/merged_refpanel.chr{chr}.benchmark.txt"
	output:
		corrs = "minimac_ref_panel/ref_panel_impacc/run{run}/merged_refpanel.chr{chr}.snp_correlations.csv", # This will contain all of the correlations for each base pair of the assay/run/chromosome
		#corrs = "imp_acc/{sample}.chr{chr}.snp_correlations.csv",
		#acc = "minimac_ref_panel/ref_panel_impacc/run{run}/merged_refpanel.chr{chr}.txt" #This file is appended to with each chromosome whose accuracy is calculated, but this can't be a valid output because it doesn't have all the wildcards in it.
	shell:
		"(python bin/vcf_impacc.py {input.true} {input.imputed} {params.acc} {output.corrs})>{log}"


rule vcf_convert:
	input:
		haps = haps_runlocator
	params:
		#inprefix = "eagle_phased_assays/run2/{sample}.chr{chr}.phased",
		inprefix = haps_runprefix,
		oprefix = "phased_assay_vcf/run{run}/{sample}.chr{chr}.phased"
	log:
		"logs/hap_leg/run{run}/{sample}.chr{chr}.phased.log"
	benchmark:
		"benchmarks/vcf_convert/run{run}/{sample}.chr{chr}.phased.benchmark.txt"
	output:
		vcf = temp("phased_assay_vcf/run{run}/{sample}.chr{chr}.phased.vcf"),
		log = "phased_assay_vcf/logs/run{run}/{sample}.chr{chr}.phased.log"
	shell:
		"(shapeit -convert --input-haps {params.inprefix} --output-log {output.log} --output-vcf {output.vcf}) > {log}"

rule run_minimac:
	input:
		ref = "minimac_ref_panel/ref_panel_chrom_vcf/run{run}/merged_refpanel.chr{chr}.phased.vcf",
		study = "phased_assay_vcf/run{run}/{sample}.chr{chr}.phased.vcf"
	params:
		oprefix = "minimac_imputed/run{run}/{sample}.chr{chr}.imputed",
		chrom = "{chr}"
	threads: 8
	log:
		"logs/minimac/{sample}.{run}.chr{chr}.phased.log"
	benchmark:
		"benchmarks/minimac/run{run}/{sample}.chr{chr}.phased.benchmark.txt"
	output:
		vcf = "minimac_imputed/run{run}/{sample}.chr{chr}.imputed.dose.vcf.gz"
	shell:
		"(Minimac3-omp --refHaps {input.ref} --haps {input.study} --myChromosome {params.chrom} --cpu 8 --prefix {params.oprefix}) > {log}" #Also consider throwing in the --hapOutput to output phased haplotype file
