SAMPLES = ['hol_testset.SNP50.788', 'hol_testset.GGPLD.788', 'hol_testset.F250.197', 'hol_testset.HD.197']

rule targ:
	input:
		targ = expand("imp_acc/{sample}.run{run}.chr{chr}.snp_correlations.csv", sample = SAMPLES, run = 1, chr = list(range(1,30)))
		#("imp_acc/{sample}.run{run}.txt:q", sample = SAMPLES, run = )

# rule impute2_vcf:
# 	input:
# 		gen = "impute2_chromosome/{sample}.run{run}.chr{chr}.phased.imputed.gen",
# 		sample = "impute2_chromosome/{sample}.chr{chr}.phased.sample"
# 	params:
# 		oprefix = "impute2_vcf/{sample}.run{run}.chr{chr}.imputed",
# 		chrom = "{chr}"
# 	log:
# 		"logs/impute2_vcf/{sample}.run{run}.chr{chr}.txt"
# 	benchmark:
# 		"benchmarks/impute2_vcf/{sample}.run{run}.chr{chr}.benchmark.txt"
# 	output:
# 		vcf = "impute2_vcf/{sample}.run{run}.chr{chr}.imputed.vcf"
# 	shell:
# 		"(plink --gen {input.gen} --sample {input.sample} --cow --keep-allele-order --oxford-single-chr {params.chrom} --recode vcf --out {params.oprefix})>{log}"
#
rule truth_chromsplit:
	input:
		gen = "accuracy_test/merged_refs/F250_HD_merged.1970.gen",
		sample = "accuracy_test/merged_refs/F250_HD_merged.1970.sample"
	params:
		oprefix = "accuracy_test/merged_refs/F250_HD_merged.chr{chr}",
		chrom = "{chr}"
	log:
		"logs/truth_chromsplit/F250_HD_merged.chr{chr}.txt"
	benchmark:
		"benchmarks/truth_chromsplit/F250_HD_merged.chr{chr}.benchmark.txt"
	output:
		vcf = "accuracy_test/merged_refs/F250_HD_merged.chr{chr}.vcf",
		frq = "accuracy_test/merged_refs/F250_HD_merged.chr{chr}.frq",
		map = "accuracy_test/merged_refs/F250_HD_merged.chr{chr}.map"
	shell:
		"(plink --gen {input.gen} --sample {input.sample} --cow --keep-allele-order --freq --chr {params.chrom} --recode vcf --out {params.oprefix}; plink --gen {input.gen} --sample {input.sample} --chr {params.chrom} --cow --keep-allele-order --recode --out {params.oprefix}; )>{log}"


rule imp_acc:
	input:
		true = "accuracy_test/merged_refs/F250_HD_merged.chr{chr}.vcf",
		imputed = "impute2_vcf/{sample}.run{run}.chr{chr}.imputed.vcf",
		#imputed = "impute2_vcf/{sample}.chr{chr}.imputed.vcf",
		frq = "accuracy_test/merged_refs/F250_HD_merged.chr{chr}.frq",
		map = "accuracy_test/merged_refs/F250_HD_merged.chr{chr}.map"
	params:
		chrom = "{chr}",
		acc = "imp_acc/{sample}.txt"
	log:
		#"logs/imp_acc/{sample}.chr{chr}.txt"
		"logs/imp_acc/{sample}.run{run}.chr{chr}.txt"
	benchmark:
		#"benchmarks/imp_acc/{sample}.chr{chr}.benchmark.txt"
		"benchmarks/imp_acc/{sample}.run{run}.chr{chr}.benchmark.txt"
	output:
		corrs = "imp_acc/{run}/{sample}.run{run}.chr{chr}.snp_correlations.csv", # This will contain all of the correlations for each base pair of the assay/run/chromosome
		#corrs = "imp_acc/{sample}.chr{chr}.snp_correlations.csv",
		acc = "imp_acc/{run}/{sample}.run{run}.txt" #This file is appended to with each chromosome whose accuracy is calculated
	shell:
		"(python bin/vcf_impacc.py {input.true} {input.imputed} {params.acc} {output.corrs})>{log}"
