from itertools import tee
def parwise(iterable):
	a,b = tee(iterable)
	next(b, None)
	return zip(a,b)

def chunks(end):
	return [str(i[0]+1)+ ' ' + str(i[1]) for i in parwise(range(0,end,5000000))]

rangedict = {'1':chunks(158322647),'2':chunks(136914030),'3':chunks(121412973), '4':chunks(120786530), '5':chunks(121186724), '6':chunks(119454666), '7':chunks(112628884), '8':chunks(113380773), '9':chunks(105701306), '10':chunks(104301732), '11':chunks(107282960), '12':chunks(91155612), '13':chunks(84230359), '14':chunks(84628243), '15':chunks(85272311), '16':chunks(81720984), '17':chunks(75149392), '18':chunks(65999195), '19':chunks(64044783), '20':chunks(71992748), '21':chunks(71594139), '22':chunks(61379134), '23':chunks(52467978), '24':chunks(62685898), '25':chunks(43879707), '26':chunks(51680135), '27':chunks(45402893), '28':chunks(46267578), '29':chunks(51504286)} #, '30':chunks(143032828)}
#describes the length of each chromsome so taht it can be chunked into equal sized chunks.

IMPGENS = ['777962.170619.411.200_B','777962.170619.2779.200_A', '227234.170619.1994.200_A']

def chrchunker(WC):
	return rangedict[WC.chr][int(WC.chunk)]

rule impute:
	input:
		targ = expand("impute4_seq_chromosome/run{run}/{sample}.chr{chr}.gen.gz", run = 26, sample = IMPGENS, chr = list(range(20,30)))

rule chip_impute2: #for parralel phasing
	input:
		hap = "/CIFS/MUG01_N/deckerje/tnr343/170519_hol/Imputation_Project/impute2_refpanel/run2/merged_refpanel.chr{chr}.{chunk}.phased.hap",
		legend = "/CIFS/MUG01_N/deckerje/tnr343/170519_hol/Imputation_Project/impute2_refpanel/run2/merged_refpanel.chr{chr}.{chunk}.phased.legend",
		knownhaps = "/CIFS/MUG01_N/deckerje/tnr343/170519_hol/Imputation_Project/eagle_phased_assays/run2/snp50.chr{chr}.phased.haps",
		maps="impute_maps/imputemap.chr{chr}.map"
	params:
		chunk= chrchunker
	log:
		"logs/chip_impute2/run{run}/{sample}.chr{chr}.{chunk}.log"
	benchmark:
		"benchmarks/chip_impute2/run{run}/{sample}.chr{chr}.{chunk}.benchmark.txt"
	output:
		imputed = "impute2_imputed/run{run}/{sample}.chr{chr}.{chunk}.gen",
		haps = "impute2_imputed/run{run}/{sample}.chr{chr}.{chunk}.gen_haps",
		summary = "impute2_imputed/run{run}/{sample}.chr{chr}.{chunk}.gen_summary"
	shell:
		"(impute2 -m {input.maps} -h {input.hap} -l {input.legend} -use_prephased_g -known_haps_g {input.knownhaps} -int {params.chunk} -phase -Ne 200 -o {output.imputed}) > {log}"

rule seqref_to_hap:
	input:
		vcfgz="/CIFS/MUG01_S/schnabelr/1kbulls/run6/beaglevcf/Chr{chr}-Beagle-TauInd-Run6.vcf.gz",
		index="/CIFS/MUG01_S/schnabelr/1kbulls/run6/beaglevcf/Chr{chr}-Beagle-TauInd-Run6.vcf.gz.tbi"
	params:
		hap = "seqref_hap/Chr{chr}-Beagle-TauInd-Run6.phased.haplotypes"
	benchmark:
		"benchmarks/seqref_to_hap/Chr{chr}-Beagle-TauInd-Run6.benchmark.txt"
	log:
		"logs/vcf_to_hap/seqref_to_hap/Chr{chr}-Beagle-TauInd-Run6.log"
	output:
		legend = "seqref_hap/Chr{chr}-Beagle-TauInd-Run6.phased.legend",
		hap = temp("seqref_hap/Chr{chr}-Beagle-TauInd-Run6.phased.haplotypes"),
		sample = "seqref_hap/Chr{chr}-Beagle-TauInd-Run6.phased.samples"
	shell:
		"(bcftools convert {input.vcfgz} --haplegendsample {params.hap},{output.legend},{output.sample}) > {log}"

rule gzip_hap:
	input:
		hap = "seqref_hap/Chr{chr}-Beagle-TauInd-Run6.phased.haplotypes"
	benchmark:
		"benchmarks/seqref_to_hap/Chr{chr}-Beagle-TauInd-Run6.benchmark.txt"
	log:
		"logs/vcf_to_hap/seqref_to_hap/Chr{chr}-Beagle-TauInd-Run6.log"
	output:
		hap = "seqref_hap/Chr{chr}-Beagle-TauInd-Run6.phased.haplotypes.gz",
	shell:
		"(gzip -c {input.hap} > {output.hap})"

rule run_impute4_seq: #for parralel phasing
	input:
		hap = "seqref_hap/Chr{chr}-Beagle-TauInd-Run6.phased.haplotypes.gz",
		legend = "seqref_hap/Chr{chr}-Beagle-TauInd-Run6.phased.legend",
		knownhaps = "impute2_imputed/run{run}/{sample}.chr{chr}.{chunk}.gen_haps",
		maps = "impute_maps/imputemap.chr{chr}.map"
	params:
		chunk= chrchunker,
		oprefix = "impute4_seq/run{run}/{sample}.chr{chr}.{chunk}"
	log:
		"logs/impute4_seq/run{run}/{sample}.chr{chr}.{chunk}.log"
	benchmark:
		"benchmarks/impute4_seq/run{run}/{sample}.chr{chr}.{chunk}.benchmark.txt"
	output:
		imputed=temp("impute4_seq/run{run}/{sample}.chr{chr}.{chunk}.gen.gz")
	shell:
		"(impute4.r265.1 -m {input.maps} -h {input.hap} -l {input.legend} -g {input.knownhaps} -int {params.chunk} -Ne 20000 -no_maf_align -o_gz -o {params.oprefix}) > {log}"

rundict = {}
for Run in range(35):
	run= str(Run)
	sampledict = {}
	for sample in IMPGENS:
		filedict = {}
		for chr in rangedict.keys():
			chunkcounter=-1
			flist = []
			for chunk in rangedict.get(chr):
				chunkcounter = chunkcounter+1
				file = 'impute4_seq/run'+ run + '/'+sample + '.chr'+chr+'.'+str(chunkcounter)+'.gen.gz' #need to edit this to accept run as a wildcard
				flist.append(file)
				filedict[chr]=flist
		sampledict[sample] = filedict
	rundict[run] = sampledict

def chrfiles(chrom):
	return rundict[chrom.run][chrom.sample][chrom.chr]

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
		cat = "impute4_seq_chromosome/run{run}/{sample}.chr{chr}.gen.gz"
	shell:
		"(cat {input.chunks} > {output.cat}) > {log}"

rule plink_convert
