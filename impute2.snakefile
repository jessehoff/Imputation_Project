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

IMPREFS = ['f250','hd']
IMPGENS = ['snp50','ggpld','f250','hd']


rule impute:
	input:
		targ = expand("impute2_chromosome/{sample}.chr{chr}.phased.imputed.gen", sample = IMPGENS, run = 2, chr = list(range(29,30)))

include: "phasing.snakefile"
include: "shapeit.snakefile"

rule impute2_refpanel:
	input:
		hap = expand("vcf_to_hap/{sample}.chr{{chr}}.phased.haplotypes", sample = IMPREFS),
		legend = expand("vcf_to_hap/{sample}.chr{{chr}}.phased.legend", sample = IMPREFS),
		maps="impute_maps/imputemap.chr{chr}.map"
	params:
		chunk = chrchunker,
		oprefix = "impute2_refpanel/merged_refpanel.chr{chr}.{chunk}.phased"
	benchmark:
		"benchmarks/impute2_refpanel/merged_refpanel.chr{chr}.{chunk}.phased.benchmark.txt"
	log:
		"logs/refpanel_impute/merged_refpanel.chr{chr}.phased.log"
	output:
		refhap = "impute2_refpanel/merged_refpanel.chr{chr}.{chunk}.phased.hap",
		refleg = "impute2_refpanel/merged_refpanel.chr{chr}.{chunk}.phased.legend"
	shell:
		"(impute2 -merge_ref_panels_output_ref {params.oprefix} -m {input.maps} -h {input.hap} -l {input.legend} -int {params.chunk} -Ne 200 -o {params.oprefix}) > {log}"

rule run_impute2_run2: #for parralel phasing
	input:
		hap = "impute2_refpanel/merged_refpanel.chr{chr}.{chunk}.phased.hap",	#How are we supposed to expand over a function? Not sure if we can make this as smart as we want it to be?
		legend = "impute2_refpanel/merged_refpanel.chr{chr}.{chunk}.phased.legend",
		knownhaps = "vcf_to_haps/{sample}.chr{chr}.phased.haps",
		maps="impute_maps/imputemap.chr{chr}.map"
	params:
		chunk= chrchunker
	log:
		"logs/impute2/{sample}.chr{chr}.{chunk}.log"
	benchmark:
		"benchmarks/impute2/{sample}.chr{chr}.{chunk}.benchmark.txt"
	output:
		imputed="impute2_imputed/{sample}.chr{chr}.{chunk}.gen",
		phasedimputed="impute2_imputed/{sample}.chr{chr}.{chunk}.phased.gen_haps"
		summary="impute2_imputed/{sample}.chr{chr}.{chunk}.phased.impute2_summary"
	shell:
		"(impute2 -merge_ref_panels -m {input.maps} -h {input.hap} -l {input.legend} -use_prephased_g -known_haps_g {input.knownhaps} -int {params.chunk} -phase -Ne 200 -o {output.imputed}) > {log}"

sampledict = {}
for sample in IMPGENS:
	filedict = {}
	for chr in rangedict.keys():
		chunkcounter=-1
		flist = []
		for chunk in rangedict.get(chr):
			chunkcounter = chunkcounter+1
			file = 'impute2_imputed/'+sample + '.chr'+chr+'.'+str(chunkcounter)+'.phased.gen_haps'
			flist.append(file)
			filedict[chr]=flist
	sampledict[sample] = filedict

def chrfiles(WC):
	return sampledict[WC.sample][WC.chr]

rule cat_chunks:
	input:
		chunks = chrfiles
	benchmark:
		"benchmarks/cat_chunks/{sample}.chr{chr}.benchmark.txt"
	log:
		"logs/cat_chunks/{sample}.chr{chr}.log"
	output:
		cat = temp("impute2_chromosome/{sample}.chr{chr}.phased.imputed.haps")
	shell:
		"(cat {input.chunks} > {output.cat}) > {log}"#; cp eagle_phased_assays/*.run2*.sample impute2_chromosome/"

#maybe needed?
 #bcftools query -l ./vcf_to_assays/hol_testset.F250.197.1.chr10.vcf > ./vcf_to_assays/hol_testset.F250.197.1.samples

#perl ./bin/vcf2impute_legend_haps.pl -vcf vcf_to_assays/hol_testset.GGPLD.788.1.chr25.vcf -leghap vcf_to_haps/hol_testset.GGPLD.788.1.chr25 -chr 25
