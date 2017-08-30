DATA =['f250', 'ggpld', 'hd', 'snp50'] #new file names -- these are files that have had ref-alt conversions
#include: 'phasing.snakefile'
rule shape_targ:
	input:
		targ = expand("shapeit_phased_assays/run{run}/{sample}.chr{chr}.phased.haps", run =4, sample = DATA, chr = list(range(1,30)))

run_dict = {'1':'merged_chrsplit', '3':'merged_chrsplit', '2':'assay_chrsplit/', '4':'assay_chrsplit/', '9':'assay_chrsplit/', '14':'assay_chrsplit/'}

def runchoice(WC):
	r = WC.run
	chrom = WC.chr
	if r == '3':
		location = run_dict[r] + '/run' + WC.run + '/hol_testset.merge.chr' + chrom +'.bed'
		#print(location)
	if r == '4':
		location = run_dict[r] + WC.sample + '.list1.chr' + chrom + '.bed'
	if r == '9':
		location = run_dict[r] + WC.sample + '.list2.chr' + chrom + '.bed'
	if r == '14':
		location = run_dict[r] + WC.sample + '.list3.chr' + chrom + '.bed'
	return location

rule run_shapeit_assay:
	input:
		bed = runchoice
	params:
		inprefix = "assay_chrsplit/{sample}.list1.chr{chr}",
		oprefix = "shapeit_phased_assays/run{run}/{sample}.chr{chr}.phased",
		map = "impute_maps/imputemap.chr{chr}.map"
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
		"(shapeit -B {params.inprefix} -M {params.map} --output-log {output.log} --thread 8 --effective-size 200 -O {params.oprefix}) > {log}"

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







def findmergelist(wc):
	listdict = {'3':'1'}
	run = wc.run
	list = listdict[run]
	locate = 'assay_chrsplit/hol_testset.list' + list+'.chr'+wc.chr+'.txt'
	#print(locate)
	return locate

rule make_shapeit_merge_list: #makes a merge list of raw genotypes, doesn't need to have a run because its just going to be a list related thing.
	input:
		filelist = expand("assay_chrsplit/{assay}.list{{list}}.chr{{chr}}.bed", assay= DATA )
	log:
		"logs/make_merge_list/list{list}.txt"
	output:
		"assay_chrsplit/hol_testset.list{list}.chr{chr}.txt"
	shell:
		"python ./bin/merge_file_maker.py {input.filelist} {output}"

rule shapeit_merged_chrsplit:
	input:
		findmergelist
	params:
		outprefix="merged_chrsplit/run{run}/hol_testset.merge.chr{chr}"
	benchmark:
		"benchmarks/merged_chrsplit/run{run}/hol_testset.merg.chr{chr}.txt"
	log:
		"logs/merged_chrsplit/run{run}/hol_testset.merge.chr{chr}.log"
	output:
		bed="merged_chrsplit/run{run}/hol_testset.merge.chr{chr}.bed",
		animalset = "merged_chrsplit/run{run}/hol_testset.merge.chr{chr}.fam"
	shell:
		"(plink --merge-list {input} --nonfounders --real-ref-alleles --cow --make-bed --out {params.outprefix}) > {log}"

rule shapeit_merged:
	input:
		bed = runchoice,
		map = "impute_maps/imputemap.chr{chr}.map"
	params:
		inprefix="merged_chrsplit/run{run}/hol_testset.merge.chr{chr}",
		oprefix="shapeit_merged/run{run}/hol_testset.merge.chr{chr}.phased"
	threads: 10
	priority: 30
	benchmark:
		"benchmarks/shapeit_merged/run{run}/hol_testset.merge.chr{chr}.benchmark.txt"
	log:
		"logs/shapeit_merged/run{run}/hol_testset.merge.chr{chr}.log"
	output:
		sample = "shapeit_merged/run{run}/hol_testset.merge.chr{chr}.phased.sample",
		haps = "shapeit_merged/run{run}/hol_testset.merge.chr{chr}.phased.haps",
		log = "shapeit_phased_assays/run{run}/hol_testset.merge.chr{chr}.phased.log"
	shell:
		"(shapeit -B {params.inprefix} -M {input.map} --output-log {output.log} --thread 8 --effective-size 200 -O {params.oprefix})> {log} "

rule shapeit_decompress:
		input:
			gzhaps = "shapeit_merged/run{run}/{sample}.chr{chr}.haps.gz"
		log:
			"logs/decompress/run{run}/{sample}.chr{chr}.log"
		benchmark:
			"benchmarks/decompress/run{run}/{sample}.chr{chr}.benchmark.txt"
		output:
			haps = temp("shapeit_merged/run{run}/{sample}.chr{chr}.haps")
		shell:
			"(gunzip -c {input.gzhaps} > {output.haps}) > {log}"

rule shapeit_merged_vcf:
	input:
		haps="shapeit_merged/run{run}/hol_testset.merge.chr{chr}.phased.haps"
	log:
			"logs/shapeit_merged_vcf/logs/run{run}/hol_testset.merge.{chr}.log"
	benchmark:
			"benchmarks/shapeit_merged_vcf/run{run}/hol_testset.merge.{chr}.benchmark.txt"
	params:
		haps="shapeit_merged/run{run}/hol_testset.merge.chr{chr}.phased"
	output:
		vcf=temp("shapeit_merged_vcf/run{run}/hol_testset.merge.chr{chr}.phased.vcf"),
		log="shapeit_merged_vcf/logs/run{run}/hol_testset.merge.chr{chr}.log"
	shell:
		"(shapeit -convert --input-haps {params.haps} --output-log {output.log} --output-vcf {output.vcf}) > {log}"

rule shapeit_bgzip_vcf:
	input:
		vcf="shapeit_merged_vcf/run{run}/hol_testset.merge.chr{chr}.phased.vcf",
	log:
		"logs/bgzip_vcf/run{run}/hol_testset.merge.chr{chr}"
	benchmark:
		"benchmarks/bgzip_vcf/run{run}/hol_testset.merge.chr{chr}"
	output:
		vcfgz=temp("shapeit_merged_vcf/run{run}/hol_testset.merge.chr{chr}.phased.vcf.gz"),
		index=temp("shapeit_merged_vcf/run{run}/hol_testset.merge.chr{chr}.phased.vcf.gz.tbi")
	shell:
		"(bgzip {input.vcf}; tabix {output.vcfgz}) > {log};"

def pickidsforlist(wc):
	listdict = {'3':'1'}
	run = wc.run
	list = listdict[run]
	return list

rule make_shapeit_vcf_extract_lists: # this rule doesn't include the list, but that is determined by the run
	input:
		animalset = "merged_chrsplit/run{run}/hol_testset.merge.chr25.fam"
	params:
		list = pickidsforlist
	benchmark:
		"benchmarks/make_vcf_extract_lists/run{run}/{assay}.benchmark.txt"
	log:
		"logs/make_vcf_extract_lists/run{run}/{assay}.benchmark.txt"
	output:
		keep_ids = "shapeit_merged_chrsplit/run{run}/phased_{assay}.keepvcf",
		keep_maps = "shapeit_merged_chrsplit/run{run}/phased_{assay}.vcfregion"
	shell:
		"python ./bin/vcfextraction_for_joint_phase.py {input.animalset} {params.list}"



rule shapeit_vcf_per_assay: #filter the vcfs on a per assay basis
	input:
		vcfgz="shapeit_merged_vcf/run{run}/hol_testset.merge.chr{chr}.phased.vcf.gz",
		index="shapeit_merged_vcf/run{run}/hol_testset.merge.chr{chr}.phased.vcf.gz.tbi",
		keep_ids = "shapeit_merged_chrsplit/run{run}/phased_{assay}.keepvcf",
		keep_maps = "shapeit_merged_chrsplit/run{run}/phased_{assay}.vcfregion"
	benchmark:
		"benchmarks/vcf_per_assay/run{run}/{assay}.chr{chr}.benchmark.txt"
	log:
		"logs/vcf_per_assay/run{run}/{assay}.chr{chr}.log"
	output:
		vcf = temp("vcf_per_assay/run{run}/{assay}.chr{chr}.vcf"),
	shell:
		"(bcftools view {input.vcfgz} -R {input.keep_maps}  -S {input.keep_ids} -o {output.vcf}) > {log}"

rule shapeit_vcf_to_hap:
	input:
		vcf = "vcf_per_assay/run{run}/{assay}.chr{chr}.vcf",
	params:
		chr = "{chr}",
		oprefix ="vcf_to_hap/run{run}/{assay}.chr{chr}"
	benchmark:
		"benchmarks/vcf_to_hap/run{run}/{assay}.chr{chr}.benchmark.txt"
	log:
		"logs/vcf_to_hap/run{run}/{assay}.chr{chr}.log"
	output:
		legend = temp("vcf_to_hap/run{run}/{assay}.chr{chr}.phased.legend"),
		haps = temp("vcf_to_hap/run{run}/{assay}.chr{chr}.phased.haplotypes"),
		sample = temp("vcf_to_hap/run{run}/{assay}.chr{chr}.phased.samples")
	shell:
		"(bcftools convert {input.vcf} --haplegendsample {output.haps},{output.legend},{output.sample}) > {log}" #Updated these to use the naming conventions that the shapeit tool outputs the run2 hap,leg,samples, so naming is consistent in Impute2.

rule shapeit_vcf_to_haps: #doesn't approrpriately name "haps" haps
	input:
		vcf = "vcf_per_assay/run{run}/{assay}.chr{chr}.vcf",
	params:
		chr = "{chr}",
		oprefix ="vcf_to_haps/run{run}/{assay}.chr{chr}.phased"
	benchmark:
		"benchmarks/vcf_to_haps/run{run}/{assay}.chr{chr}.benchmark.txt"
	log:
		"logs/vcf_to_haps/run{run}/{assay}.chr{chr}.log"
	output:
		hap = temp("vcf_to_haps/run{run}/{assay}.chr{chr}.phased.haps"),
		sample = temp("vcf_to_haps/run{run}/{assay}.chr{chr}.phased.sample")
	shell:
		"(bcftools convert {input.vcf} --hapsample {output.hap},{output.sample} ) > {log}"
