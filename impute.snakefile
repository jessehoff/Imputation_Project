from itertools import tee
def parwise(iterable):
        a,b = tee(iterable)
        next(b, None)
        return zip(a,b)

def chunks(end):
        return [str(i[0]+1)+ ' ' + str(i[1]) for i in parwise(range(0,end,5000000))]

rangedict = {'1':chunks(158322647),'2':chunks(136914030),'3':chunks(121412973), '4':chunks(120786530), '5':chunks(121186724), '6':chunks(119454666), '7':chunks(112628884), '8':chunks(113380773), '9':chunks(105701306), '10':chunks(104301732), '11':chunks(107282960), '12':chunks(91155612), '13':chunks(84230359), '14':chunks(84628243), '15':chunks(85272311), '16':chunks(81720984), '17':chunks(75149392), '18':chunks(65999195), '19':chunks(64044783), '20':chunks(71992748), '21':chunks(71594139), '22':chunks(61379134), '23':chunks(52467978), '24':chunks(62685898), '25':chunks(43879707), '26':chunks(51680135), '27':chunks(45402893), '28':chunks(46267578), '29':chunks(51504286)}

def chrchunker(WC):
        return rangedict[WC.chr][int(WC.chunk)]

DATA =['227234.170112.325.GGPF250', '26504.170112.3126.GGPLDV3', '58336.170112.315.SNP50A', '58336.170112.335.SNP50B', '58336.170112.3399.SNP50C', '76999.170112.3498.GGP90KT','777962.170127.483.HD']
IMPREFS1 =['58336.170112.315.SNP50A', '58336.170112.335.SNP50B', '58336.170112.3399.SNP50C', '76999.170112.3498.GGP90KT']
IMPGENS1 =['26504.170112.3126.GGPLDV3']
IMPREFS2 =['227234.170112.325.GGPF250', '777962.170127.483.HD']
IMPGENS2 =['227234.170112.325.GGPF250', '58336.170112.315.SNP50A', '58336.170112.335.SNP50B', '58336.170112.3399.SNP50C', '76999.170112.3498.GGP90KT','777962.170127.483.HD']

targfiles=[]
for chr in rangedict.keys():
    chunkcounter=0
    for chunk in rangedict.get(chr):
        chunkcounter = chunkcounter
        file = 'imprun1/imp_round1.chr'+chr+'.'+str(chunkcounter)+'.phased.impute2'
        targfiles.append(file)
#  "imprun1/imp_round1.chr{chr}.{chunk}.phased.impute2"
rule targ:
	input:
        	targ=targfiles[:1]
        #targ = expand("impute_input/{sample}.chr{chr}.phased.gen.gz", sample = DATA, chr = list(range(1,30)))


rule decompress:
	input:
		gzhaps = "eagle_phased_assays/{sample}.chr{chr}.phased.haps.gz"
	log:
		"snake_logs/decompress/{sample}.chr{chr}.phased.log"
	output:
		haps = "eagle_phased_assays/{sample}.chr{chr}.phased.haps"
	shell:
		"(gunzip -c {input.gzhaps} > {output.haps})>{log}" #rm eagle_phased_assays/*.haps.gz) > {log}"

rule convert:
	input:
		haps= "eagle_phased_assays/{sample}.chr{chr}.phased.haps",
		sample = "eagle_phased_assays/{sample}.chr{chr}.phased.sample"
	params:
		inprefix = "eagle_phased_assays/{sample}.chr{chr}.phased",
		oprefix = "impute_vcf/{sample}.chr{chr}.phased"
	log:
		"snake_logs/convert/{sample}.chr{chr}.phased.log"
	benchmark:
		"filter_benchmarks/convert/{sample}.chr{chr}.phased.benchmark.txt"
	output:
		vcf = "impute_vcf/{sample}.chr{chr}.phased.vcf"
	shell:
		"(shapeit -convert --input-haps {params.inprefix} --output-vcf {output.vcf}) > {log}"

#Takes VCF input for each assay, and creates corresponding haps, legend, sample, and gen files for each that will become inputs for impute2
rule impute_input:
    input:
        vcf = "impute_vcf/{sample}.chr{chr}.phased.vcf"
    params:
        oprefix = "impute_input/{sample}.chr{chr}.phased",
        genfix = "impute_input/{sample}.chr{chr}.phased.gen"
    log:
        "snake_logs/impute_input_creation/{sample}.chr{chr}.log"
    benchmark:
        "filter_benchmarks/impute_input_creation/{sample}.chr{chr}.phased.benchmark.txt"
    output:
        legend = "impute_input/{sample}.chr{chr}.phased.legend.gz",
        hap = "impute_input/{sample}.chr{chr}.phased.hap.gz",
        gen = "impute_input/{sample}.chr{chr}.phased.gen.gz",
        #sample = "impute_input/{sample}.chr{chr}.phased.gen.sample"
    shell:
        "(perl bin/vcf2impute_legend_haps.pl -vcf {input.vcf} -leghap {params.oprefix} -chr {chr}; perl bin/vcf2impute_gen.pl -vcf {input.vcf} -gen {params.genfix} -chr {chr}) > {log}"
#def gener(WC):
#	genslist = []
#	for i in IMPGENS1:
#		genslist.append('impute_input/' + i + '.chr' + WC.chr + '.phased.gen.gz')
#	return genslist

rule run_impute2_round1:
        input:
                haps=expand("impute_input/{sample}.chr{{chr}}.phased.hap.gz", sample=IMPREFS1) ,
                legend=expand("impute_input/{sample}.chr{{chr}}.phased.legend.gz", sample=IMPREFS1),
                gens= expand("impute_input/{sample}.chr{{chr}}.phased.gen.gz", sample = IMPGENS1),
                maps="impute_maps/imputemap.chr{chr}.map"
        params:
                chunk= chrchunker
        log:
                "logs/imprun1/imp_round1.chr{chr}.{chunk}.log"
        benchmark:
                "benchmarks/imprun1/imp_round1.chr{chr}.{chunk}.log"
        output:
                imputed="imprun1/imp_round1.chr{chr}.{chunk}.phased.impute2",
                summary="imprun1/imp_round1.chr{chr}.{chunk}.phased.impute2_summary"
        shell:
                "(impute2 -merge_ref_panels -m {input.maps} -h {input.haps} -l {input.legend} -use_prephased_g -known_haps_g {input.gens} -int {params.chunk} -Ne 100 -o {output.imputed}) > {log}"
