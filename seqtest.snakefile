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

rule target:
	input:
		targ = expand("vcf_to_haps/{sample}.chr{chr}.phased.haps", sample = ["ggpld", "snp50", "f250", "hd"], chr = list(range(1,30)))
rule merged_chrsplit:
	input:
		mergelist = "dataprepper/F250_HD_merged.1970.bed"
	params:
		outprefix="merged_chrsplit/F250_HD_merged.1970.chr{chr}",
		inprefix = "dataprepper/F250_HD_merged.1970",
		chrom = "{chr}"
	benchmark:
		"benchmarks/chrsplit/F250_HD_merged.1970.chr{chr}.txt"
	log:
		"logs/chrsplit/F250_HD_merged.1970.chr{chr}.txt"
	output:
		bed="merged_chrsplit/F250_HD_merged.1970.chr{chr}.bed",
	shell:
		"(plink --bfile {params.inprefix} --nonfounders --real-ref-alleles --cow --chr {params.chrom} --make-bed --out {params.outprefix}) > {log}"

rule eagle_phasing:
	input:
		bed="merged_chrsplit/F250_HD_merged.1970.chr{chr}.bed"
	params:
		inprefix="merged_chrsplit/F250_HD_merged.1970.chr{chr}",
		oprefix="eagle_merged/F250_HD_merged.1970.chr{chr}"
	threads: 10
	priority: 30
	benchmark:
		"benchmarks/eagle_merged/F250_HD_merged.1970.chr{chr}.benchmark.txt"
	log:
		"logs/eagle_merged/F250_HD_merged.1970.chr{chr}.log"
	output:
		sample = "eagle_merged/F250_HD_merged.1970.chr{chr}.sample",
		haps = "eagle_merged/F250_HD_merged.1970.chr{chr}.haps.gz"
	shell:
		"(eagle --bfile={params.inprefix}  --geneticMapFile=USE_BIM --maxMissingPerSnp .99  --maxMissingPerIndiv .99 --numThreads 10 --outPrefix {params.oprefix})> {log}"

rule decompress:
	input:
		gzhaps = "eagle_merged/F250_HD_merged.1970.chr{chr}.haps.gz"
	log:
		"logs/decompress/eagle_merged/F250_HD_merged.1970.chr{chr}.log"
	benchmark:
		"benchmarks/decompress/eagle_merged/F250_HD_merged.1970.chr{chr}.benchmark.txt"
	output:
		haps = temp("eagle_merged/F250_HD_merged.1970.chr{chr}.haps")
	shell:
		"(gunzip -c {input.gzhaps} > {output.haps}) > {log}"

rule eagle_merged_vcf:
	input:
		haps="eagle_merged/F250_HD_merged.1970.chr{chr}.haps"
	log:
			"logs/eagle_merged_vcf/logs/F250_HD_merged.1970.chr{chr}.log"
	benchmark:
			"benchmarks/eagle_merged_vcf/F250_HD_merged.1970.chr{chr}.benchmark.txt"
	params:
		haps="eagle_merged/F250_HD_merged.1970.chr{chr}"
	output:
		vcf=temp("eagle_merged_vcf/F250_HD_merged.1970.chr{chr}.phased.vcf"),
		log="eagle_merged_vcf/logs/F250_HD_merged.1970.chr{chr}.log"
	shell:
		"(shapeit -convert --input-haps {params.haps} --output-log {output.log} --output-vcf {output.vcf}) > {log}"

rule bgzip_vcf:
	input:
		vcf="eagle_merged_vcf/F250_HD_merged.1970.chr{chr}.phased.vcf",
	log:
		"logs/bgzip_vcf/F250_HD_merged.1970.chr{chr}.log"
	benchmark:
		"benchmarks/bgzip_vcf/F250_HD_merged.1970.chr{chr}.benchmark.txt"
	output:
		vcfgz=temp("eagle_merged_vcf/F250_HD_merged.1970.chr{chr}.phased.vcf.gz"),
		index=temp("eagle_merged_vcf/F250_HD_merged.1970.chr{chr}.phased.vcf.gz.tbi")
	shell:
		"(bgzip {input.vcf}; tabix {output.vcfgz}) > {log}"

rule vcf_per_assay: #filter the vcfs on a per assay basis
	input:
		vcfgz="eagle_merged_vcf/F250_HD_merged.1970.chr{chr}.phased.vcf.gz",
		index="eagle_merged_vcf/F250_HD_merged.1970.chr{chr}.phased.vcf.gz.tbi",
		keep_maps= "dataprepper/{sample}_snps_umd3.1.txt",
		keep_ids = "dataprepper/brd_1kbulls_list.txt"
	benchmark:
		"benchmarks/vcf_per_assay/{sample}.chr{chr}.benchmark.txt"
	log:
		"logs/vcf_per_assay/{sample}.chr{chr}.log"
	output:
		vcf = temp("vcf_per_assay/{sample}.chr{chr}.vcf")
	shell:
		"(bcftools view {input.vcfgz} -R {input.keep_maps}  -S {input.keep_ids} -o {output.vcf}) > {log}"

rule vcf_to_haps: #doesn't approrpriately name "haps" haps
	input:
		vcf = "vcf_per_assay/{sample}.chr{chr}.vcf",
	params:
		chr = "{chr}",
		oprefix ="vcf_to_haps/{sample}.chr{chr}.phased"
	benchmark:
		"benchmarks/vcf_to_haps/{sample}.chr{chr}.benchmark.txt"
	log:
		"logs/vcf_to_haps/{sample}.chr{chr}.log"
	output:
		hap = "vcf_to_haps/{sample}.chr{chr}.phased.haps",
		sample = "vcf_to_haps/{sample}.chr{chr}.phased.sample"
	shell:
		"(bcftools convert {input.vcf} --hapsample {output.hap},{output.sample} ) > {log}"

rule subset_run6:
	input:
		vcfgz="/CIFS/MUG01_S/schnabelr/1kbulls/run6/beaglevcf/Chr{chr}-Beagle-TauInd-Run6.vcf.gz",
		index="/CIFS/MUG01_S/schnabelr/1kbulls/run6/beaglevcf/Chr{chr}-Beagle-TauInd-Run6.vcf.gz.tbi",
		keep_ids = "dataprepper/brd_1kbulls_list.txt"
	benchmark:
		"benchmarks/seqref/seqref.chr{chr}.benchmark.txt"
	log:
		"logs/seqref/seqref.chr{chr}.log"
	output:
		vcf = "seqref/seqref_nobrd.chr{chr}.vcf.gz",
		tabix = "seqref/seqref_nobrd.chr{chr}.vcf.gz.tbi"
	shell:
		"(bcftools view {input.vcfgz} -S {input.keep_ids} -O b -o {output.vcf}; tabix {output.vcf}) > {log}"

rule vcf_to_hap:
	input:
		vcf = "seqref/seqref_nobrd.chr{chr}.vcf.gz",
	params:
		chr = "{chr}",
		oprefix ="seqref/seqref_nobrd.chr{chr}",
		hapsout = "seqref/seqref_nobrd.chr{chr}.phased.haplotypes"
	benchmark:
		"benchmarks/vcf_to_hap/seqref_nobrd.chr{chr}.benchmark.txt"
	log:
		"logs/vcf_to_hap/seqref_nobrd.chr{chr}.log"
	output:
		legend = "seqref/seqref_nobrd.chr{chr}.phased.legend",
		haps = "seqref/seqref_nobrd.chr{chr}.phased.haplotypes.gz",
		sample = "seqref/seqref_nobrd.chr{chr}.phased.samples"
	shell:
		"(bcftools convert {input.vcf} --haplegendsample {params.hapsout},{output.legend},{output.sample}; pigz {params.hapsout}) > {log}"

rule run_impute4_seq: #for parralel phasing
	input:
		legend = "seqref/seqref_nobrd.chr{chr}.phased.legend",
		hap = "seqref/seqref_nobrd.chr{chr}.phased.haplotypes.gz",
		sample = "seqref/seqref_nobrd.chr{chr}.phased.samples",
		knownhaps = "vcf_to_haps/{sample}.chr{chr}.phased.haps",
		maps = "impute_maps/imputemap.chr{chr}.map"
	params:
		chunk= chrchunker,
		oprefix = "impute4_seq/{sample}.chr{chr}.{chunk}"
	log:
		"logs/impute4_seq/{sample}.chr{chr}.{chunk}.log"
	benchmark:
		"benchmarks/impute4_seq/{sample}.chr{chr}.{chunk}.benchmark.txt"
	output:
		imputed="impute4_seq/{sample}.chr{chr}.{chunk}.gen.gz"
	shell:
		"(impute4.r265.1 -m {input.maps} -h {input.hap} -l {input.legend} -g {input.knownhaps} -int {params.chunk} -Ne 20000 -no_maf_align -o_gz -o {params.oprefix}) > {log}"

rundict = {}
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
