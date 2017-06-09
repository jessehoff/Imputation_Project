#SAMPLES = ['hol_testset.SNP50.788', 'hol_testset.GGPLD.788', 'hol_testset.F250.197', 'hol_testset.HD.197']
SAMPLES = ['f250','snp50','ggpld','hd']

rule targ:
	input:
		#targ = expand("imp_acc/{run}/{sample}.run{run}.chr{chr}.snp_correlations.csv", sample = SAMPLES, run = [1,2], chr = list(range(1,30)))
		#("imp_acc/{sample}.run{run}.txt:q", sample = SAMPLES, run = )
		targ = expand("imp_acc_visualization/run{run}/{sample}.run{run}.chr{chr}.histogram.png", sample = SAMPLES, run=[1,2], chr = list(range(1,30)))

rule impute2_vcf:
	input:
		gen = "impute2_chromosome/{sample}.run{run}.chr{chr}.phased.imputed.gen",
		#sample = "impute2_chromosome/{sample}.chr{chr}.run{run}.phased.sample"
		sample = "impute2_chromosome/{sample}.run{run}.chr{chr}.sample"
	params:
		oprefix = "impute2_vcf/{sample}.run{run}.chr{chr}.imputed",
		chrom = "{chr}"
	log:
		"logs/impute2_vcf/{sample}.run{run}.chr{chr}.txt"
	benchmark:
		"benchmarks/impute2_vcf/{sample}.run{run}.chr{chr}.benchmark.txt"
	output:
		vcf = "impute2_vcf/{sample}.run{run}.chr{chr}.imputed.vcf"
	shell:
		"(plink --gen {input.gen} --sample {input.sample} --cow --real-ref-alleles --oxford-single-chr {params.chrom} --recode vcf --out {params.oprefix})>{log}"

rule merge_ref: #makes a merge list of raw genotypes
	input:
		HD = "correct_sex/777962.170519.1970.HD.bed",
		F250 = "correct_sex/227234.170519.1970.GGPF250.bed"
	params:
		bfile = "--bfile correct_sex/227234.170519.1970.GGPF250",
		bmerge = "--bmerge correct_sex/777962.170519.1970.HD",
		ofile = "--out merge_ref/hol_testset.F250_HD_merged.1970"
	output:
		"merge_ref/hol_testset.F250_HD_merged.1970.bed"
	shell:
		"plink  {params.bfile} {params.bmerge}  --cow --nonfounders --real-ref-alleles {params.ofile} --make-bed"


rule truth_chromsplit:
	input:
		bed = "merge_ref/hol_testset.F250_HD_merged.1970.bed",
	params:
		oprefix = "ref_vcfs/F250_HD_merged.chr{chr}",
		chrom = "{chr}",
		iprefix = "merge_ref/hol_testset.F250_HD_merged.1970"
	log:
		"logs/truth_chromsplit/F250_HD_merged.chr{chr}.txt"
	benchmark:
		"benchmarks/truth_chromsplit/F250_HD_merged.chr{chr}.benchmark.txt"
	output:
		vcf = "ref_vcfs/F250_HD_merged.chr{chr}.vcf",
		frq = "ref_vcfs/F250_HD_merged.chr{chr}.frq",
		map = "ref_vcfs/F250_HD_merged.chr{chr}.map",
		ped = temp("ref_vcfs/F250_HD_merged.chr{chr}.ped")
	shell:
		"(plink --bfile {params.iprefix}  --cow --real-ref-alleles --freq --chr {params.chrom} --recode vcf --out {params.oprefix}; plink --bfile {params.iprefix} --cow --real-ref-alleles --recode --out {params.oprefix})>{log}"


rule ref_to_df:
	input:
		true = "ref_vcfs/F250_HD_merged.chr{chr}.vcf",
	output:
		truepickle = "ref_vcfs/F250_HD_merged.chr{chr}.pickle"
	log:
		"logs/ref_to_df/F250_HD{chr}.txt"
	benchmark:
		"benchmarks/ref_to_df/F250_HD{chr}.txt"
	shell:
		"(python bin/vcf_ref_todf.py {input.true}  {output.truepickle})>{log}"


rule imp_acc:
	input:
		true = "ref_vcfs/F250_HD_merged.chr{chr}.pickle",
		imputed = "impute2_vcf/{sample}.run{run}.chr{chr}.imputed.vcf",
		#imputed = "impute2_vcf/{sample}.chr{chr}.imputed.vcf",
		frq = "ref_vcfs/F250_HD_merged.chr{chr}.frq",
	params:
		chrom = "{chr}",
		acc = "imp_acc/{run}/{sample}.txt"
	log:
		"logs/imp_acc/{sample}.run{run}.chr{chr}.txt"
	benchmark:
		#"benchmarks/imp_acc/{sample}.chr{chr}.benchmark.txt"
		"benchmarks/imp_acc/{sample}.run{run}.chr{chr}.benchmark.txt"
	output:
		corrs = "imp_acc/{run}/{sample}.run{run}.chr{chr}.snp_correlations.csv", # This will contain all of the correlations for each base pair of the assay/run/chromosome
		#corrs = "imp_acc/{sample}.chr{chr}.snp_correlations.csv",
		#acc = "imp_acc/{run}/{sample}.run{run}.txt" #This file is appended to with each chromosome whose accuracy is calculated, but this can't be a valid output because it doesn't have all the wildcards in it.
	shell:
		"(python bin/vcf_impacc.py {input.true} {input.imputed} {params.acc} {output.corrs})>{log}"

rule imp_acc_visualization:
	input:
		corrs = "imp_acc/{run}/{sample}.run{run}.chr{chr}.snp_correlations.csv",
		frq = "ref_vcfs/F250_HD_merged.chr{chr}.frq",
		map = "ref_vcfs/F250_HD_merged.chr{chr}.map"
	params:
		acc = "imp_acc_visualization/{run}/{sample}.txt"
	log:
		"logs/imp_acc_visualization/{sample}.run{run}.chr{chr}.txt"
	benchmark:
		"benchmarks/imp_acc_visualization/{sample}.run{run}.chr{chr}.benchmark.txt"
	output:
		hist = "imp_acc_visualization/run{run}/{sample}.run{run}.chr{chr}.histogram.png",
		scatter = "imp_acc_visualization/run{run}/{sample}.run{run}.chr{chr}.scatter.png",
		line = "imp_acc_visualization/run{run}/{sample}.run{run}.chr{chr}.line.png",
		combo = "imp_acc_visualization/run{run}/{sample}.run{run}.chr{chr}.combo.png"
	shell:
		"(python bin/impacc_visualization.py {input.corrs} {input.frq} {input.map} {output.hist} {output.scatter} {output.line} {output.combo}) > {log}"
