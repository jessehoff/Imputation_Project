SAMPLES = ['hol_testset.SNP50.788', 'hol_testset.GGPLD.788', 'hol_testset.F250.197', 'hol_testset.HD.197']

rule targ:
	input:
		targ = expand("imp_acc/{sample}.chr{chr}.txt", sample = SAMPLES, chr = list(range(1,30)))


rule impute2_vcf:
	input:
		gen = "impute2_chromosome/{sample}.chr{chr}.phased.imputed.gen",
		sample = "impute2_chromosome/{sample}.chr{chr}.phased.sample"
	params:
		oprefix = "impute2_vcf/{sample}.chr{chr}.imputed",
		chrom = "{chr}"
	log:
		"logs/impute2_vcf/{sample}.chr{chr}.txt"
	benchmark:
		"benchmarks/impute2_vcf/{sample}.chr{chr}.benchmark.txt"
	output:
		vcf = "impute2_vcf/{sample}.chr{chr}.imputed.vcf"
	shell:
		"(plink --gen {input.gen} --sample {input.sample} --cow --keep-allele-order --oxford-single-chr {params.chrom} --recode vcf --out {params.oprefix})>{log}"

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
		vcf = "accuracy_test/merged_refs/F250_HD_merged.chr{chr}.vcf"
	shell:
		"(plink --gen {input.gen} --sample {input.sample} --cow --keep-allele-order --chr {params.chrom} --recode vcf --out {params.oprefix})>{log}"

rule imp_acc:
	input:
		true = "accuracy_test/merged_refs/F250_HD_merged.chr{chr}.vcf",
		imputed = "impute2_vcf/{sample}.chr{chr}.imputed.vcf"
	log:
		"logs/imp_acc/{sample}.chr{chr}.txt"
	benchmark:
		"benchmarks/imp_acc/{sample}.chr{chr}.benchmark.txt"
	output:
		acc = "imp_acc/{sample}.chr{chr}.txt"
	shell:
		"(python bin/vcf_impacc.py {input.true} {input.imputed} {output.acc})>{log}"
