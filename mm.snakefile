IMPGENS =['hol_testset.SNP50.788', 'hol_testset.GGPLD.788', 'hol_testset.F250.197', 'hol_testset.HD.197']
IMPREFS =['hol_testset.F250.197', 'hol_testset.HD.197']
F250 = ['hol_testset.F250.197']
HD = ['hol_testset.HD.197']

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

rule hap_leg:
	input:
		haps = "eagle_phased_assays/{sample}.run{run}.chr{chr}.phased.haps",
		sample = "eagle_phased_assays/{sample}.run{run}.chr{chr}.phased.sample"
	params:
		inprefix = "eagle_phased_assays/{sample}.run{run}.chr{chr}.phased",
		oprefix = "reference_panel_accuracy/impute_input/{sample}.run{run}.chr{chr}.phased"
	# log:
	# 	"snake_logs/hap_leg/{sample}.run{run}.chr{chr}.phased.log"
	# benchmark:
	# 	"filter_benchmarks/hap_leg/{sample}.run{run}.chr{chr}.phased.benchmark.txt"
	output:
		hap = "reference_panel_accuracy/impute_input/{sample}.run{run}.chr{chr}.phased.haplotypes",
		leg = "reference_panel_accuracy/impute_input/{sample}.run{run}.chr{chr}.phased.legend",
		log = "reference_panel_accuracy/impute_input/logs/{sample}.run{run}.chr{chr}.phased.log"
	shell:
		"(shapeit -convert --input-haps {params.inprefix} --output-log {output.log} --output-ref {params.oprefix}) > {log}"

rule impute2_refpanel:
	input:
		hap = expand("reference_panel_accuracy/impute_input/{sample}.run{{run}}.chr{{chr}}.phased.haplotypes", sample = IMPREFS),
		leg = expand("reference_panel_accuracy/impute_input/{sample}.run{{run}}.chr{{chr}}.phased.legend", sample = IMPREFS),
		haps = "eagle_phased_assays/{sample}.run{run}.chr{chr}.phased.haps",
		maps="impute_maps/imputemap.chr{chr}.map"
	params:
		chunk = chrchunker,
		oprefix = "reference_panel_accuracy/ref_panels/logs/merged_refpanel.run{run}.chr{chr}"
	output:
		refhap = "reference_panel_accuracy/ref_panels/merged_refpanel.run{run}.chr{chr}.phased.hap",
		refleg = "reference_panel_accuracy/ref_panels/merged_refpanel.run{run}.chr{chr}.phased.legend"
	shell:
		"(impute2 -merge_ref_panels_output_ref -m {input.maps} -h {input.hap} -l {input.legend} -use_prephased_g -known_haps_g {input.knownhaps} -int {params.chunk} -Ne 200 -o {params.oprefix}) > {log}"


rule vcf_convert:
	input:
		haps = "eagle_phased_assays/{sample}.chr{chr}.phased.haps",
		sample = "eagle_phased_assays/{sample}.chr{chr}.phased.sample"
	params:
		inprefix = "eagle_phased_assays/{sample}.chr{chr}.phased",
		oprefix = "phased_assay_vcf/{sample}.chr{chr}.phased"
	log:
		"snake_logs/hap_leg/{sample}.chr{chr}.phased.log"
	benchmark:
		"filter_benchmarks/vcf_convert/{sample}.chr{chr}.phased.benchmark.txt"
	output:
		vcf = "phased_assay_vcf/{sample}.chr{chr}.phased.vcf",
		log = "phased_assay_vcf/logs/{sample}.chr{chr}.phased.log"
	shell:
		"(shapeit -convert --input-haps {params.inprefix} --output-log {output.log} --output-vcf {params.oprefix}) > {log}"

rule ref_up:
	input:
		hd = expand("phased_assay_vcf/{sample}.chr{chr}.phased.vcf", sample = HD)
		f250 = expand("phased_assay_vcf/{sample}.chr{chr}.phased.vcf", sample = F250)
	params:
		oprefix = "mm_reference/{sample}.chr{chr}.mm.ref"
	log:
		"snake_logs/f250_to_hd/{sample}.chr{chr}.phased.log"
	benchmark:
		"filter_benchmarks/f250_to_hd/{sample}.chr{chr}.phased.benchmark.txt"
	output:
		vcf = "mm_reference/{sample}.chr{chr}.dose.vcf"
