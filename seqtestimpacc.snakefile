#include: "seqtest.snakefile"
DATA = ['merged']
rule target:
	input:
		targ = expand("imp_acc_variable/problem_samples_removed/{sample}.chr{chr}.variable_snps.txt", sample = DATA, chr = list(range(1,30)))
	shell:
		"python bin/cleanup.py"
rule decompress_imp:
	input:
		gen = "impute4_seq_chromosome/{sample}.chr{chr}.gen.gz",
	#	sample = "eagle_merged/F250_HD_merged.1970.chr{chr}.sample"
	output:
		gen = temp("impute4_seq_chromosome/{sample}.chr{chr}.gen")
	shell:
		"gunzip {input.gen}; sed -i 's/Chr//g' {output.gen}"

rule impute4_seq_vcf:
	input:
		gen = "impute4_seq_chromosome/{sample}.chr{chr}.gen",
		sample = "vcf_to_haps/{sample}.chr{chr}.phased.sample"
	params:
		oprefix = "impute4_seq_vcf/{sample}.chr{chr}.imputed",
		chrom = "{chr}"
	log:
		"logs/impute4_seq_vcf/{sample}.chr{chr}.txt"
	benchmark:
		"benchmarks/impute4_seq_vcf/{sample}.chr{chr}.benchmark.txt"
	output:
		gen = "impute4_seq_chromosome/{sample}.chr{chr}.gen.gz",
		vcf = temp("impute4_seq_vcf/{sample}.chr{chr}.imputed.vcf")
	shell:
		"(bcftools convert -G {input.gen},{input.sample} -o {output.vcf};pigz {input.gen})>{log}"

rule ref_pop_frqs:
	input:
		vcfgz="/CIFS/MUG01_S/schnabelr/1kbulls/run6/beaglevcf/Chr{chr}-Beagle-TauInd-Run6.vcf.gz"
	params:
		oprefix = "impacc/seqref.chr{chr}"
	log:
		"logs/ref_pop_frqs/seqref_truth.txt"
	benchmark:
		"benchmarks/ref_pop_frqs/seqref_truth.benchmark.txt"
	output:
		frq = "impacc/seqref.chr{chr}.frq"
	shell: "(plink --vcf {input.vcfgz}  --cow --real-ref-alleles --freq --out {params.oprefix})>{log}"


rule decompres_ref:
	input:
		true = "impacc/seqref_truth.chr{chr}.vcf.gz"
	params:
		vcf = "impacc/seqref_truth.chr{chr}.vcf"
	log:
		"logs/ref_to_df/seqref_truth.chr{chr}.txt"
	benchmark:
		"benchmarks/ref_to_df/seqref_truth.chr{chr}.txt"
	output:
		vcf = "impacc/seqref_truth.chr{chr}.vcf"
	shell:
		"(gunzip {input.true})>{log}"

rule imp_acc:
	input:
		true = "impacc/seqref_truth.chr{chr}.vcf",
		imputed = "impute4_seq_vcf/{sample}.chr{chr}.imputed.vcf"
	log:
		"logs/imp_acc/{sample}.chr{chr}.txt"
	benchmark:
		"benchmarks/imp_acc/{sample}.chr{chr}.benchmark.txt"
	output:
		corrs = "imp_acc/problem_samples_removed/{sample}.chr{chr}.snp_correlations.csv", # This will contain all of the correlations for each base pair of the assay/run/chromosome_truth.chr{chr}.vcf.gz"
	shell:
		"(python bin/sequence_impacc.py {input.true} {input.imputed} {output.corrs}; pigz {input.true})>{log}"

rule imp_acc_variable:
	input:
		true = "impacc/seqref_truth.chr{chr}.vcf",
		imputed = "impute4_seq_vcf/{sample}.chr{chr}.imputed.vcf"
	log:
		"logs/imp_acc/{sample}.chr{chr}.txt"
	benchmark:
		"benchmarks/imp_acc/{sample}.chr{chr}.benchmark.txt"
	output:
		corrs = "imp_acc_variable/problem_samples_removed/{sample}.chr{chr}.variable_snps.txt", # This will contain all of the correlations for each base pair of the assay/run/chromosome_truth.chr{chr}.vcf.gz"
	shell:
		"(python bin/sequence_variable.py {input.true} {input.imputed} {output.corrs}; pigz {input.true})>{log}"
