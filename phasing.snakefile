DATA =['f250', 'ggpld', 'hd', 'snp50'] #new file names -- these are files that have had ref-alt conversions
hd_or_f250  = {'snp50':"correct_sex/777962.170519.1970.HD",'f250':"correct_sex/227234.170519.1970.GGPF250", 'ggpld':"correct_sex/777962.170519.1970.HD", 'hd':"correct_sex/777962.170519.1970.HD"}

rule targ:
	input:
		hd_or_f250  = {'snp50':"correct_sex/777962.170519.1970.HD",'f250':"correct_sex/227234.170519.1970.GGPF250", 'ggpld':"correct_sex/777962.170519.1970.HD", 'hd':"correct_sex/777962.170519.1970.HD"}

rule assay_chrsplit:
	input:
		bed = expand("downsample/{sample}.bed", assay = DATA),
	params:
		inprefix = "downsample/{sample}",
		oprefix = "assay_chrsplit/{sample}.chr{chr}",
		chr = "{chr}"
	benchmark:
		"benchmarks/assay_chrsplit/{sample}.chr{chr}.txt"
	log:
		"logs/eagle_split_chromosomes/{sample}.chr{chr}.log"
	output:
		bed = "assay_chrsplit/{sample}.chr{chr}.bed",
		bim = "assay_chrsplit/{sample}.chr{chr}.bim",
		fam = "assay_chrsplit/{sample}.chr{chr}.fam",
		log = "assay_chrsplit/{sample}.chr{chr}.log"
	shell:
		"(plink --bfile {params.inprefix}  --real-ref-alleles --chr {params.chr} --make-bed  --nonfounders --cow --out {params.oprefix})> {log}"

rule eagle_phased_assays:
	input:
		bed = "assay_chrsplit/{sample}.chr{chr}.bed"
	params:
		inprefix = "assay_chrsplit/{sample}.chr{chr}",
		oprefix = "eagle_phased_assays/{sample}.chr{chr}.phased"
	benchmark:
		"benchmarks/eagle_phased_assays/{sample}.chr{chr}.benchmark.txt"
	log:
		"logs/eagle_phased_assays/{sample}.chr{chr}.log"
	threads: 8
	priority: 100
	output:
		sample = "eagle_phased_assays/{sample}.chr{chr}.phased.sample",
		haps = "eagle_phased_assays/{sample}.chr{chr}.phased.haps.gz"
	shell:
		"(eagle --bfile={params.inprefix} --geneticMapFile=USE_BIM --maxMissingPerSnp 1 --maxMissingPerIndiv 1 --numThreads 8 --outPrefix {params.oprefix}) > {log}"

rule decompress_single_chrom:
	input:
		gzhaps = "eagle_phased_assays/{sample}.chr{chr}.phased.haps.gz"
	benchmark:
		"benchmarks/decompress/{sample}.chr{chr}.benchmark.txt"
	log:
		"logs/decompress/{sample}.chr{chr}.log"
	output:
		haps = temp("eagle_phased_assays/{sample}.chr{chr}.phased.haps")
	shell:
		"(gunzip -c {input.gzhaps} > {output.haps}) > {log}"

rule hap_leg:
	input:
		haps = "eagle_phased_assays/{sample}.chr{chr}.phased.haps",
		sample = "eagle_phased_assays/{sample}.chr{chr}.phased.sample"
	params:
		inprefix = "eagle_phased_assays/{sample}.chr{chr}.phased",
		oprefix = "impute_input/{sample}.chr{chr}.phased"
	log:
		"logs/hap_leg/{sample}.chr{chr}.phased.log"
	benchmark:
		"benchmarks/hap_leg/{sample}.chr{chr}.phased.benchmark.txt"
	output:
		hap = temp("impute_input/{sample}.chr{chr}.phased.haplotypes"),
		leg = temp("impute_input/{sample}.chr{chr}.phased.legend"),
		log = "impute_input/logs/{sample}.chr{chr}.phased.log"
	shell:
		"(shapeit -convert --input-haps {params.inprefix} --output-log {output.log} --output-ref {params.oprefix}) > {log}"

rule make_merge_list: #How are we doing this on a repeated basis? 
	input:
		filelist = expand("assay_chrsplit/{assay}.chr{{chr}}.bed", assay= DATA )
	log:
		"logs/make_merge_list/list{list}.txt"
	output:
		"assay_chrsplit/all_assays.chr{chr}.txt"
	shell:
		"python ./bin/merge_file_maker.py {input.filelist} {output}"

rule merged_chrsplit:
	input:
		filelist = "assay_chrsplit/all_assays.chr{chr}.txt"
	params:
		oprefix="merged_chrsplit/merged_assays.chr{chr}"
	benchmark:
		"benchmarks/merged_chrsplit/merged_assays.chr{chr}.txt"
	log:
		"logs/merged_chrsplit/merged_assays.chr{chr}.log"
	output:
		bed="merged_chrsplit/merged_assays.chr{chr}.bed",
		animalset = "merged_chrsplit/merged_assays.chr{chr}.fam"
	shell:
		"(plink --merge-list {input} --nonfounders --real-ref-alleles --cow --make-bed --out {params.oprefix}) > {log}"

rule eagle_merged:
	input:
		bed = "merged_chrsplit/merged_assays.chr{chr}.bed"
	params:
		bed="merged_chrsplit/merged_assays.chr{chr}",
		out="eagle_merged/merged_assays.chr{chr}"
	threads: 10
	priority: 30
	benchmark:
		"benchmarks/eagle_merged/merged_assays.chr{chr}.benchmark.txt"
	log:
		"logs/eagle_merged/merged_assays.chr{chr}.log"
	output:
		sample = "eagle_merged/merged_assays.chr{chr}.phased.sample",
		haps = "eagle_merged/merged_assays.chr{chr}.phased.haps.gz"
	shell:
		"(eagle --bfile={params.bed}  --geneticMapFile=USE_BIM --maxMissingPerSnp .99  --maxMissingPerIndiv .99 --numThreads 10 --outPrefix {params.out})> {log} "

rule decompress: #Is this still needed?
	input:
		gzhaps = "eagle_merged/merged_assays.chr{chr}.haps.gz"
	benchmark:
		"benchmarks/decompress/merged_assays.chr{chr}.benchmark.txt"
	log:
		"logs/decompress/merged_assays.chr{chr}.log"
	output:
		haps = temp("eagle_merged/merged_assays.chr{chr}.haps")
	shell:
		"(gunzip -c {input.gzhaps} > {output.haps}) > {log}"

rule eagle_merged_vcf:
	input:
		haps="eagle_merged/merged_assays.chr{chr}.haps"
	benchmark:
		"benchmarks/eagle_merged_vcf/merged_assays.{chr}.benchmark.txt"
	log:
		"logs/eagle_merged_vcf/logs/merged_assays.{chr}.log"
	params:
		haps="eagle_merged/merged_assays.chr{chr}"
	output:
		vcf=temp("eagle_merged_vcf/merged_assays.chr{chr}.phased.vcf"),
		log="eagle_merged_vcf/logs/merged_assays.chr{chr}.log"
	shell:
		"(shapeit -convert --input-haps {params.haps} --output-log {output.log} --output-vcf {output.vcf}) > {log}"

rule bgzip_vcf:
	input:
		vcf="eagle_merged_vcf/merged_assays.chr{chr}.phased.vcf",
	log:
		"logs/bgzip_vcf/merged_assays.chr{chr}"
	benchmark:
		"benchmarks/bgzip_vcf/merged_assays.chr{chr}"
	output:
		vcfgz=temp("eagle_merged_vcf/merged_assays.chr{chr}.phased.vcf.gz"),
		index=temp("eagle_merged_vcf/merged_assays.chr{chr}.phased.vcf.gz.tbi")
	shell:
		"(bgzip {input.vcf}; tabix {output.vcfgz}) > {log};"


rule make_vcf_extract_lists: # This rule should be much easier when now that we're not trying to make lists talk back and forth.  Directly reference
	input:
		bim = "assay_chrsplit/{sample}.chr{chr}.bim",
		fam = "assay_chrsplit/{sample}.chr{chr}.fam"
	benchmark:
		"benchmarks/make_vcf_extract_lists/{sample}.chr{chr}.benchmark.txt"
	log:
		"logs/make_vcf_extract_lists/{sample}.chr{chr}.log"
	output:
		keep_ids = "merged_chrsplit/extract_lists/{sample}.chr{chr}.keepvcf",
		keep_snps = "merged_chrsplit/extract_lists/{sample}.chr{chr}.vcfregion"
	shell:
		"python vcf_extraction_maker.py {input.bim} {input.fam} {output.keep_snps} {output.keep_ids}"


rule vcf_per_assay: #filter the vcfs on a per chromosome per assay basis
	input:
		vcfgz="eagle_merged_vcf/merged_assays.chr{chr}.phased.vcf.gz",
		index="eagle_merged_vcf/merged_assays.chr{chr}.phased.vcf.gz.tbi",
		keep_ids = "merged_chrsplit/extract_lists/{sample}.chr{chr}.keepvcf",
		keep_snps = "merged_chrsplit/extract_lists/{sample}.chr{chr}.vcfregion"
	benchmark:
		"benchmarks/vcf_per_assay/{sample}.chr{chr}.benchmark.txt"
	log:
		"logs/vcf_per_assay/{sample}.chr{chr}.log"
	output:
		vcf = temp("vcf_per_assay/{sample}.chr{chr}.vcf"),
	shell:
		"(bcftools view {input.vcfgz} -R {input.keep_snps} -S {input.keep_ids} -o {output.vcf}) > {log}"

rule vcf_to_hap:
	input:
		vcf = "vcf_per_assay/{sample}.chr{chr}.vcf",
	params:
		chr = "{chr}",
		oprefix ="vcf_to_hap/{sample}.chr{chr}"
	benchmark:
		"benchmarks/vcf_to_hap/{sample}.chr{chr}.benchmark.txt"
	log:
		"logs/vcf_to_hap/{sample}.chr{chr}.log"
	output:
		legend = temp("vcf_to_hap/{sample}.chr{chr}.phased.legend"),
		hap = temp("vcf_to_hap/{sample}.chr{chr}.phased.haplotypes"),
		sample = temp("vcf_to_hap/{sample}.chr{chr}.phased.samples")
	shell:
		"(bcftools convert {input.vcf} --haplegendsample {output.hap},{output.legend},{output.sample}) > {log}" #Uses the naming conventions that the shapeit tool outputs for hap,leg,samples

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
		hap = temp("vcf_to_haps/{sample}.chr{chr}.phased.haps"),
		sample = temp("vcf_to_haps/{sample}.chr{chr}.phased.sample")
	shell:
		"(bcftools convert {input.vcf} --hapsample {output.hap},{output.sample} ) > {log}"
