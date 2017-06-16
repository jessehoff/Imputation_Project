#SAMPLES = ['hol_testset.SNP50.788', 'hol_testset.GGPLD.788', 'hol_testset.F250.197', 'hol_testset.HD.197']
SAMPLES = ['f250','snp50','ggpld','hd']

rule impacc:
	input:
		#targ = expand("imp_acc/{run}/{sample}.run{run}.chr{chr}.snp_correlations.csv", sample = SAMPLES, run = [1,2], chr = list(range(1,30)))
		#("imp_acc/{sample}.run{run}.txt:q", sample = SAMPLES, run = )
		#targ = expand("imp_acc_visualization/run{run}/{sample}.run{run}.chr{chr}.histogram.png", sample = SAMPLES, run=[1,2], chr = list(range(1,30)))
		#targ = expand("ref_vcfs/F250_HD_merged.chr{chr}.pickle", chr = list(range(1,30)))
		#targ = expand("minimac_imp_acc/{run}/{sample}.run{run}.chr{chr}.snp_correlations.csv", run = 2, sample = SAMPLES, chr = list(range(1,30)))
		targ = expand("imp_acc/run{run}/{sample}.mafcorr.csv", run = 1, sample = SAMPLES)
		#targ = expand("imp_acc/run{run}/{sample}.lowmafcorr.png", run = 6, sample = SAMPLES)
#include: "mm.snakefile"
include: "impute2.snakefile"

#rule sample_file_move:
#	input:
#		samp = "/run{run}/{sample}.chr{chr}.phased.sample"	#Need to create some sort of a sample file mover function, also naming convention of phased/imputed data needs to be carried through
#		samp = "/run
#	params:
#		oprefix = "impute2_chromosome/run{run}/"
#	log:
#		"logs/sample_file_move/run{run}/{sample}.chr{chr}.log"
#	output:
#		samp = "impute2_chromosome/run{run}/{sample}.chr{chr}.phased.sample"
#	shell:
#		"cp {input.samp} {params.oprefix}"

#sample = "vcf_to_hap/run{run}/{assay}.chr{chr}.phased.samples"

def samplefinder(WC):
	rundict = {'6':"vcf_to_haps",'1':"vcf_to_haps"}
	r = WC.run
	chrom = WC.chr
	if rundict[r] == "vcf_to_haps":
		location = rundict[r] + '/run' + r+'/' + WC.sample+'.chr' + chrom+'.phased.sample'
	return location

rule impute2_vcf:
	input:
		gen = "impute2_chromosome/run{run}/{sample}.chr{chr}.phased.imputed.gen",
		sample = samplefinder
	params:
		oprefix = "impute2_vcf/run{run}/{sample}.chr{chr}.imputed",
		chrom = "{chr}"
	log:
		"logs/impute2_vcf/run{run}/{sample}.chr{chr}.txt"
	benchmark:
		"benchmarks/impute2_vcf/run{run}/{sample}.chr{chr}.benchmark.txt"
	output:
		vcf = "impute2_vcf/run{run}/{sample}.chr{chr}.imputed.vcf"
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
		#idslist = "--keep dataprepper/{assay}_ids.list{list}.txt"
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
		imputed = "impute2_vcf/run{run}/{sample}.chr{chr}.imputed.vcf",
		#imputed = "impute2_vcf/{sample}.chr{chr}.imputed.vcf",
		frq = "ref_vcfs/F250_HD_merged.chr{chr}.frq",
	params:
		chrom = "{chr}",
		acc = "imp_acc/run{run}/{sample}.txt"
	log:
		"logs/imp_acc/run{run}/{sample}.chr{chr}.txt"
	benchmark:
		#"benchmarks/imp_acc/{sample}.chr{chr}.benchmark.txt"
		"benchmarks/imp_acc/run{run}/{sample}.chr{chr}.benchmark.txt"
	output:
		corrs = "imp_acc/run{run}/{sample}.chr{chr}.snp_correlations.csv", # This will contain all of the correlations for each base pair of the assay/run/chromosome
		#corrs = "imp_acc/{sample}.chr{chr}.snp_correlations.csv",
		#acc = "imp_acc/run{run}/{sample}.run{run}.txt" #This file is appended to with each chromosome whose accuracy is calculated, but this can't be a valid output because it doesn't have all the wildcards in it.
	shell:
		"(python bin/vcf_impacc.py {input.true} {input.imputed} {params.acc} {output.corrs})>{log}"

rule imp_acc_visualization:
	input:
		corrs = "imp_acc/run{run}/{sample}.chr{chr}.snp_correlations.csv",
		frq = "ref_vcfs/F250_HD_merged.chr{chr}.frq",
		map = "ref_vcfs/F250_HD_merged.chr{chr}.map"
	params:
		acc = "imp_acc/run{run}/visualization/{sample}.txt"
	log:
		"logs/imp_acc_visualization/run{run}/{sample}.chr{chr}.txt"
	benchmark:
		"benchmarks/imp_acc_visualization/run{run}/{sample}.chr{chr}.benchmark.txt"
	output:
		hist = "imp_acc/run{run}/visualization/{sample}.chr{chr}.histogram.png",
		scatter = "imp_acc/run{run}/visualization/{sample}.chr{chr}.scatter.png",
		line = "imp_acc/run{run}/visualization/{sample}.chr{chr}.line.png",
		combo = "imp_acc/run{run}/visualization/{sample}.chr{chr}.combo.png"
	shell:
		"(python bin/impacc_visualization.py {input.corrs} {input.frq} {input.map} {output.hist} {output.scatter} {output.line} {output.combo}) > {log}"

rule all_chrom_impacc:
	input:
		corrs = expand("imp_acc/run{{run}}/{{sample}}.chr{chr}.snp_correlations.csv", chr = list(range(1,30)))
	params:
		corrprefix = "imp_acc/run{run}/{sample}.chr",
		frq = "accuracy_test/merged_refs/F250_HD_merged.1970.frq",
		mapfile = "accuracy_test/merged_refs/F250_HD_merged.1970.map",
		master = "imp_acc/master_impacc.txt"
	log:
		"logs/all_chrom_impacc/run{run}/{sample}.txt"
	benchmark:
		"benchmarks/all_chrom_impacc/run{run}/{sample}.benchmark.txt"
	output:
		corrout = "imp_acc/run{run}/{sample}.mafcorr.csv",
		fig = "imp_acc/run{run}/{sample}.mafcorr.png",
		lowmaf_fig = "imp_acc/run{run}/{sample}.lowmafcorr.png"
	shell:
		"(python bin/allchrom_impacc.py {params.corrprefix} {params.frq} {params.mapfile} {output.corrout} {output.fig} {output.lowmaf_fig} {params.master})"
