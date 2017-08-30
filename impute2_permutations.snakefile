rule permutations:
	input:
		targ = expand("imp_acc/run{run}/visualization/{sample}.chr{chr}.lowmaf.png", run = list(range(32,81)), sample = 'snp50', chr = 20)

IMPREFS = ['f250','hd']
IMPGENS = ['snp50','ggpld','f250','hd']

nedict = {'32':'50', '33':'100', '34':'250', '35':'500', '36':'1000', '37':'2500', '38':'5000', '39':'10000', '40':'20000', '41':'25000', '42':'29999'}
def nefinder(WC):
	return nedict[WC.run]

callthreshdict = {'43':'0.1', '44':'0.2', '45':'0.3', '46':'0.4', '47':'0.5', '48':'0.6', '49':'0.7', '50':'0.8', '51':'0.9', '52':'1.0'}
def callthreshfinder(WC):
	return callthreshdict[WC.run]

khapsdict = {'53':'50', '54':'100', '55':'200', '56':'300', '57':'400', '58':'500', '59':'600', '60':'700', '61':'800'}
def khapsfinder(WC):
	return khapsdict[WC.run]

kdict = {'62':'20', '63':'50', '64':'80', '65':'100', '66':'200', '67':'300', '68':'400','69':'800', '70':'1000'}
def kfinder(WC):
	return kdict[WC.run]

iterdict = {'71':'11', '72':'20', '73':'50', '74':'100', '75':'150'}
def iterfinder(WC):
	return iterdict[WC.run]

burndict = {'76':'0', '77':'5', '78':'10', '79':'20', '80':'29'}
def burnfinder(WC):
	return burndict[WC.run]

# chunkrundict = {'81':'100000', '82':'500000', '83':'1000000', '84':'2000000', '85':'3000000', '86':'4000000', '87':'5000000', '88':'1000000'}
# def chunksize(WC):
# 	return chunkrundict[WC.run]


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
		hap = expand("vcf_to_hap/run1/{sample}.chr{{chr}}.phased.haplotypes", sample = IMPREFS),
		legend = expand("vcf_to_hap/run1/{sample}.chr{{chr}}.phased.legend", sample = IMPREFS),
		#maps = "impute_maps/recombination.chr{chr}.map"
		maps="impute_maps/imputemap.chr{chr}.map"
	log:
		"logs/refpanel_impute/run{run}/merged_refpanel.chr{chr}.phased.log"
	params:
		chunk = chrchunker,
		oprefix = "impute2_refpanel/run{run}/merged_refpanel.chr{chr}.{chunk}.phased",
		#burn = burnfinder
		#iter = iterfinder
		#k = kfinder
		#khaps = khapsfinder
		#callthresh = callthreshfinder
		#ne = nefinder
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
		knownhaps = "vcf_to_haps/run1/{sample}.chr{chr}.phased.haps",
		#maps = "impute_maps/recombination.chr{chr}.map"
		maps="impute_maps/imputemap.chr{chr}.map"
	params:
		chunk= chrchunker,
		#burn = burnfinder
		#iter = iterfinder
		#k = kfinder
		#khaps = khapsfinder
		#callthresh = callthreshfinder
		#ne = nefinder
	log:
		"logs/impute2/run{run}/{sample}.chr{chr}.{chunk}.log"
	benchmark:
		"benchmarks/impute2/run{run}/{sample}.chr{chr}.{chunk}.benchmark.txt"
	output:
		imputed="impute2_imputed/run{run}/{sample}.chr{chr}.{chunk}.impute2",
		summary="impute2_imputed/run{run}/{sample}.chr{chr}.{chunk}.impute2_summary"
	shell:
		"(impute2 -merge_ref_panels -m {input.maps} -h {input.hap} -l {input.legend} -use_prephased_g -known_haps_g {input.knownhaps} -int {params.chunk} -phase -Ne 200 -o {output.imputed}) > {log}"

# rule run_impute2_run2: #for parralel phasing
# 	input:
# 		hap = "vcf_to_hap/run1/hd.chr{chr}.phased.haplotypes",	#How are we supposed to expand over a function? Not sure if we can make this as smart as we want it to be?
# 		legend = "vcf_to_hap/run1/hd.chr{chr}.phased.legend",
# 		knownhaps = "vcf_to_haps/run1/{sample}.chr{chr}.phased.haps",
# 		#maps = "impute_maps/recombination.chr{chr}.map"
# 		maps="impute_maps/imputemap.chr{chr}.map"
# 	params:
# 		chunk= chrchunker,
# 		#burn = burnfinder
# 		#iter = iterfinder
# 		#k = kfinder
# 		#khaps = khapsfinder
# 		#callthresh = callthreshfinder
# 		#ne = nefinder
# 	log:
# 		"logs/impute2/run{run}/{sample}.chr{chr}.{chunk}.log"
# 	benchmark:
# 		"benchmarks/impute2/run{run}/{sample}.chr{chr}.{chunk}.benchmark.txt"
# 	output:
# 		imputed="impute2_imputed/run{run}/{sample}.chr{chr}.{chunk}.impute2",
# 		summary="impute2_imputed/run{run}/{sample}.chr{chr}.{chunk}.impute2_summary"
# 	shell:
# 		"(impute2 -merge_ref_panels -m {input.maps} -h {input.hap} -l {input.legend} -use_prephased_g -known_haps_g {input.knownhaps} -int {params.chunk} -fill_holes -Ne 200 --o {output.imputed}) > {log}"


rundict = {}
for Run in range(100):
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
		sample = "vcf_to_haps/run1/{sample}.chr{chr}.phased.sample"
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
