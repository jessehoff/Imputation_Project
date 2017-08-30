DATA =['f250', 'ggpld', 'hd', 'snp50']

rule target:
	input:
		targ = expand("shapeit_phased_assays/run{run}/{sample}.chr{chr}.phased.haps", run = 19, sample = DATA, chr = list(range(1,30)))

run_dict = {'17':'assay_chrsplit/', '18':'assay_chrsplit/','19':'assay_chrsplit/'}

def runchoice(WC):
	r = WC.run
	chrom = WC.chr
	if r == '17':
		location = run_dict[r] + WC.sample + '.list1.chr' + chrom + '.bed'
	if r == '18':
		location = run_dict[r] + WC.sample + '.list2.chr' + chrom + '.bed'
	if r == '19':
		location = run_dict[r] + WC.sample + '.list3.chr' + chrom + '.bed'
	return location

rule run_shapeit_assay:
	input:
		bed = runchoice,
		map = "impute_maps/imputemap.chr{chr}.map"
	params:
		inprefix = "assay_chrsplit/{sample}.list1.chr{chr}",
		oprefix = "shapeit_phased_assays/run{run}/{sample}.chr{chr}.phased"
	benchmark:
		"benchmarks/shapeit_phased_assays/run{run}/{sample}.chr{chr}.benchmark.txt"
	threads: 8
	log:
		"logs/shapeit_phased_assays/run{run}/{sample}.chr{chr}.log"
	output:
		sample = "shapeit_phased_assays/run{run}/{sample}.chr{chr}.phased.sample",
		haps = "shapeit_phased_assays/run{run}/{sample}.chr{chr}.phased.haps",
		log = "shapeit_phased_assays/run{run}/{sample}.chr{chr}.phased.log"
	shell:
		"(shapeit3.r884.1 -B {params.inprefix} -M {input.map} --output-log {output.log} --thread 8 --effective-size 200 -O {params.oprefix}) > {log}"

rule shapeit_hap_leg:
	input:
		haps = "shapeit_phased_assays/run{run}/{sample}.chr{chr}.phased.haps",
		sample = "shapeit_phased_assays/run{run}/{sample}.chr{chr}.phased.sample"
	params:
		inprefix = "shapeit_phased_assays/run{run}/{sample}.chr{chr}.phased",
		oprefix = "shapeit_phased_assays/impute_input/run{run}/{sample}.chr{chr}.phased"
	log:
		"logs/hap_leg/run{run}/{sample}.chr{chr}.phased.log"
	benchmark:
		"benchmarks/hap_leg/run{run}/{sample}.chr{chr}.phased.benchmark.txt"
	output:
		hap = "shapeit_phased_assays/impute_input/run{run}/{sample}.chr{chr}.phased.haplotypes",
		leg = "shapeit_phased_assays/impute_input/run{run}/{sample}.chr{chr}.phased.legend",
		log = "shapeit_phased_assays/impute_input/run{run}/logs/{sample}.chr{chr}.phased.log"
	shell:
		"(shapeit -convert --input-haps {params.inprefix} --output-log {output.log} --output-ref {params.oprefix}) > {log}"
