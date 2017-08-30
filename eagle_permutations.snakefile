rule permute:
	input:
		targ = expand("imp_acc/run{run}/visualization/{sample}.chr{chr}.lowmaf.png", run = list(range(101, 131)), sample = 'snp50', chr = 20)

paramdict = {'101':'--Kpbwt 100', '102':'--Kpbwt 250', '103':'--Kpbwt 500', '104':'--Kpbwt 1000', '105':'--Kpbwt 2500', '106':'--Kpbwt 5000', '107':'--Kpbwt 10000', '108':'--Kpbwt 15000', '109':'--Kpbwt 20000', '110':'--Kpbwt 30000', '111':'--pbwtIters 0', '112':'--pbwtIters 1', '113':'--pbwtIters 2', '114':'--pbwtIters 3', '115':'--pbwtIters 3','116':'--expectIBDcM  0.5', '117':'--expectIBDcM  1', '118':'--expectIBDcM  2', '119':'--expectIBDcM  3', '120':'--expectIBDcM  5', '121':'--histFactor 0', '122':'--histFactor 1', '123':'--histFactor 2', '124':'--histFactor 3', '125':'--histFactor 5', '126':'--genoErrProb 0.001', '127':'--genoErrProb 0.003', '128':'--genoErrProb 0.005', '129':'--genoErrProb 0.01', '130':'--genoErrProb 0.1'}
def paramfinder(WC):
	return paramdict[WC.run]

rule eagle_phased_assays:
	input:
		bed = "assay_chrsplit/{sample}.list1.chr{chr}.bed"
	params:
		inprefix = "assay_chrsplit/{sample}.list1.chr{chr}",
		oprefix = "eagle_phased_assays/run{run}/{sample}.chr{chr}.phased",
		perm = paramfinder
	benchmark:
		"benchmarks/eagle_phased_assays/run{run}/{sample}.chr{chr}.benchmark.txt"
	threads: 8
	priority: 100
	log:
		"logs/eagle_phased_assays/run{run}/{sample}.chr{chr}.log"
	output:
		sample = "eagle_phased_assays/run{run}/{sample}.chr{chr}.phased.sample",
		haps = "eagle_phased_assays/run{run}/{sample}.chr{chr}.phased.haps.gz"
	shell:
		"(eagle --bfile={params.inprefix} --geneticMapFile=USE_BIM --maxMissingPerSnp 1 --maxMissingPerIndiv 1 --numThreads 8 {params.perm} --outPrefix {params.oprefix}) > {log}"

rule decompress_single_chrom:
		input:
			gzhaps = "eagle_phased_assays/run{run}/{sample}.chr{chr}.phased.haps.gz"		# inprefix = listchoice,
		log:
			"logs/decompress/run{run}/{sample}.chr{chr}.log"
		benchmark:
			"benchmarks/decompress/run{run}/{sample}.chr{chr}.benchmark.txt"
		output:
			haps = temp("eagle_phased_assays/run{run}/{sample}.chr{chr}.phased.haps")
		shell:
			"(gunzip -c {input.gzhaps} > {output.haps}) > {log}"

rule hap_leg:
	input:
		haps = "eagle_phased_assays/run{run}/{sample}.chr{chr}.phased.haps",
		sample = "eagle_phased_assays/run{run}/{sample}.chr{chr}.phased.sample"
	params:
		inprefix = "eagle_phased_assays/run{run}/{sample}.chr{chr}.phased",
		oprefix = "impute_input/run{run}/{sample}.chr{chr}.phased"
	log:
		"logs/hap_leg/run{run}/{sample}.chr{chr}.phased.log"
	benchmark:
		"benchmarks/hap_leg/run{run}/{sample}.chr{chr}.phased.benchmark.txt"
	output:
		hap = temp("impute_input/run{run}/{sample}.chr{chr}.phased.haplotypes"),
		leg = temp("impute_input/run{run}/{sample}.chr{chr}.phased.legend"),
		log = "impute_input/run{run}/logs/{sample}.chr{chr}.phased.log"
	shell:
		"(shapeit -convert --input-haps {params.inprefix} --output-log {output.log} --output-ref {params.oprefix}) > {log}"

IMPREFS = ['f250','hd']
IMPGENS = ['snp50']

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

rule impute2_refpanel:
	input:
		hap = expand("impute_input/run{{run}}/{sample}.chr{{chr}}.phased.haplotypes", sample = IMPREFS),
		legend = expand("impute_input/run{{run}}/{sample}.chr{{chr}}.phased.legend", sample = IMPREFS),
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

rule run_impute2_run2:
	input:
		hap = "impute2_refpanel/run{run}/merged_refpanel.chr{chr}.{chunk}.phased.hap",	#How are we supposed to expand over a function? Not sure if we can make this as smart as we want it to be?
		legend = "impute2_refpanel/run{run}/merged_refpanel.chr{chr}.{chunk}.phased.legend",
		knownhaps = "eagle_phased_assays/run{run}/{sample}.chr{chr}.phased.haps",
		maps="impute_maps/imputemap.chr{chr}.map"
	params:
		chunk= chrchunker
	log:
		"logs/impute2/run{run}/{sample}.chr{chr}.{chunk}.log"
	benchmark:
		"benchmarks/impute2/run{run}/{sample}.chr{chr}.{chunk}.benchmark.txt"
	output:
		imputed="impute2_imputed/run{run}/{sample}.chr{chr}.{chunk}.impute2",
		summary="impute2_imputed/run{run}/{sample}.chr{chr}.{chunk}.impute2_summary"
	shell:
		"(impute2 -merge_ref_panels -m {input.maps} -h {input.hap} -l {input.legend} -use_prephased_g -known_haps_g {input.knownhaps} -int {params.chunk} -phase -Ne 200 -o {output.imputed}) > {log}"

rundict = {}
for Run in range(150):
	run= str(Run)
	sampledict = {}
	for sample in IMPGENS:
		filedict = {}
		for chr in rangedict.keys():
			chunkcounter=-1
			flist = []
			for chunk in rangedict.get(chr):
				chunkcounter = chunkcounter+1
				file = 'impute2_imputed/run'+ run + '/'+sample + '.chr'+chr+'.'+str(chunkcounter)+'.impute2' #need to edit this to accept run as a wildcard
				flist.append(file)
				filedict[chr]=flist
		sampledict[sample] = filedict
	rundict[run] = sampledict

def chrfiles(chrom):
	return rundict[chrom.run][chrom.sample][chrom.chr]


rule cat_chunks:
	input:
		chunks = chrfiles
	log:
		"logs/cat_chunks/run{run}/{sample}.chr{chr}.log"
	benchmark:
		"benchmarks/cat_chunks/run{run}/{sample}.chr{chr}.benchmark.txt"
	output:
		cat = temp("impute2_chromosome/run{run}/{sample}.chr{chr}.imputed.gen")
	shell:
		"(cat {input.chunks} > {output.cat}) > {log}"

rule impute2_vcf:
	input:
		gen = "impute2_chromosome/run{run}/{sample}.chr{chr}.imputed.gen",
		sample = "eagle_phased_assays/run{run}/{sample}.chr{chr}.phased.sample"
	params:
		oprefix = "impute2_vcf/run{run}/{sample}.chr{chr}.imputed",
		chrom = "{chr}"
	log:
		"logs/impute2_vcf/run{run}/{sample}.chr{chr}.txt"
	benchmark:
		"benchmarks/impute2_vcf/run{run}/{sample}.chr{chr}.benchmark.txt"
	output:
		vcf = temp("impute2_vcf/run{run}/{sample}.chr{chr}.imputed.vcf")
	shell:
		"(plink --gen {input.gen} --sample {input.sample} --cow --real-ref-alleles --oxford-single-chr {params.chrom} --recode vcf --out {params.oprefix})>{log}"

rule imp_acc:
	input:
		true = "ref_vcfs/F250_HD_merged.chr{chr}.pickle",
		imputed = "impute2_vcf/run{run}/{sample}.chr{chr}.imputed.vcf",
		frq = "ref_vcfs/F250_HD_merged.chr{chr}.frq",
	params:
		chrom = "{chr}",
		acc = "imp_acc/permutations.chr20.accuracies.txt"
	log:
		"logs/imp_acc/run{run}/{sample}.chr{chr}.txt"
	benchmark:
		"benchmarks/imp_acc/run{run}/{sample}.chr{chr}.benchmark.txt"
	output:
		corrs = "imp_acc/run{run}/{sample}.chr{chr}.snp_correlations.csv"
	shell:
		"(python bin/vcf_impacc.py {input.true} {input.imputed} {params.acc} {output.corrs})>{log}"

rule imp_acc_visualization:
	input:
		corrs = "imp_acc/run{run}/{sample}.chr{chr}.snp_correlations.csv",
		frq = "ref_vcfs/F250_HD_merged.chr{chr}.frq",
		map = "ref_vcfs/F250_HD_merged.chr{chr}.map"
	params:
		acc = "imp_acc/run{run}/visualization/{sample}.txt"
	log:
		"logs/imp_acc_visualization/run{run}/{sample}.chr{chr}.txt"
	benchmark:
		"benchmarks/imp_acc_visualization/run{run}/{sample}.chr{chr}.benchmark.txt"
	output:
		hist = "imp_acc/run{run}/visualization/{sample}.chr{chr}.histogram.png",
		scatter = "imp_acc/run{run}/visualization/{sample}.chr{chr}.scatter.png",
		line = "imp_acc/run{run}/visualization/{sample}.chr{chr}.line.png",
		combo = "imp_acc/run{run}/visualization/{sample}.chr{chr}.combo.png",
		lowmaf = "imp_acc/run{run}/visualization/{sample}.chr{chr}.lowmaf.png"
	shell:
		"(python bin/impacc_visualization.py {input.corrs} {input.frq} {input.map} {output.hist} {output.scatter} {output.line} {output.combo} {output.lowmaf}) > {log}"
