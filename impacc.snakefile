#SAMPLES = ['hol_testset.SNP50.788', 'hol_testset.GGPLD.788', 'hol_testset.F250.197', 'hol_testset.HD.197']
SAMPLES = ['f250','snp50','ggpld','hd']

rule impacc:
	input:
		#targ = expand("imp_acc/{run}/{sample}.run{run}.chr{chr}.snp_correlations.csv", sample = SAMPLES, run = [1,2], chr = list(range(1,30)))
		#("imp_acc/{sample}.run{run}.txt:q", sample = SAMPLES, run = )
		#targ = expand("imp_acc_visualization/run{run}/{sample}.run{run}.chr{chr}.histogram.png", sample = SAMPLES, run=[1,2], chr = list(range(1,30)))
		#targ = expand("ref_vcfs/F250_HD_merged.chr{chr}.pickle", chr = list(range(1,30)))
		#targ = expand("minimac_imp_acc/{run}/{sample}.run{run}.chr{chr}.snp_correlations.csv", run = 2, sample = SAMPLES, chr = list(range(1,30)))
		#targ = expand("imp_acc/run{run}/{sample}.mafcorr.csv", run = 5, sample = SAMPLES)
		#targ = expand("imp_acc/run{run}/visualization/{sample}.chr{chr}.combo.png", run = 10, sample = SAMPLES, chr = 28)
		targ = expand("imp_acc/run{run}/{sample}.mafcorr.csv", run = 11, sample = SAMPLES)
		#targ = expand("imp_acc/run{run}/{sample}.lowmafcorr.png", run = 6, sample = SAMPLES)
include: "mm.snakefile"
include: "impute2.snakefile"


def samplefinder(WC):
	rundict = {'6':"vcf_to_haps",'1':"vcf_to_haps",'12':"vcf_to_haps", '2':'eagle_phased_assays','13':'eagle_phased_assays', '4':"shapeit_phased_assays", '7':'eagle_phased_assays', '9':'shapeit_phased_assays'}
	r = WC.run
	chrom = WC.chr
	location = rundict[r] + '/run' + r+'/' + WC.sample+'.chr' + chrom+'.phased.sample'
	return location

def impaccscript(WC):
	scriptdict = {'1':'bin/vcf_impacc.py', '2':'bin/vcf_impacc.py', '3':'bin/vcf_impacc.py', '4':'bin/vcf_impacc.py', '5':'bin/minimac_impacc.py', '6':'bin/vcf_impacc.py', '7':'bin/vcf_impacc.py', '8':'bin/minimac_impacc.py', '9':'bin/vcf_impacc.py', '10':'bin/minimac_impacc.py', '11':'bin/minimac_impacc.py', '12':'bin/vcf_impacc.py', '13':'bin/vcf_impacc.py', '14':'bin/vcf_impacc.py', '15':'bin/minimac_impacc.py', '16':'bin/minimac_impacc.py'}
	script = scriptdict[WC.run]
	return script
# def impute2vcffinder(WC):
# 	dirdict = {'1':'impute2_vcf', '2':'impute2_vcf', '3':'impute2_vcf', '4':'impute2_vcf','6':'impute2_vcf', '7':'impute2_vcf', '9':'impute2_vcf','12':'impute2_vcf','9':'impute2_vcf'}
# 	suffdict = {'1':'.imputed.vcf','12':'.imputed.vcf', '2':'.imputed.vcf', '3':'.imputed.vcf', '4':'.imputed.vcf', '6':'.imputed.vcf', '7':'.imputed.vcf', '9':'.imputed.vcf'}
# 	location = 'impute2_vcf' + '/run' + WC.run + '/' + WC.sample + '.chr' + WC.chr + '.imputed.vcf'
# 	return location
#
def vcffinder(WC):
	dirdict = {'1':'impute2_vcf', '2':'impute2_vcf', '3':'impute2_vcf', '4':'impute2_vcf','5':'minimac_imputed', '6':'impute2_vcf', '7':'impute2_vcf', '8':'minimac_imputed', '9':'impute2_vcf', '10':'minimac_imputed', '11':'minimac_imputed', '12':'impute2_vcf','13':'impute2_vcf','14':'impute2_vcf', '15':'minimac_imputed', '16':'minimac_imputed'}
	suffdict = {'1':'.imputed.vcf', '2':'.imputed.vcf', '3':'.imputed.vcf', '4':'.imputed.vcf', '6':'.imputed.vcf', '7':'.imputed.vcf', '9':'.imputed.vcf', '5':'.imputed.dose.vcf', '8':'.imputed.dose.vcf', '10':'.imputed.dose.vcf', '11':'.imputed.dose.vcf','12':'.imputed.vcf','13':'.imputed.vcf','14':'.imputed.vcf', '15':'.imputed.dose.vcf', '16':'.imputed.dose.vcf'}
	location = dirdict[WC.run] + '/run' + WC.run + '/' + WC.sample + '.chr' + WC.chr + suffdict[WC.run]
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

rule ref_pop_frqs:
	input:
		bed = "merge_ref/hol_testset.F250_HD_merged.1970.bed",
	params:
		oprefix = "ref_vcfs/F250_HD_merged",
		iprefix = "merge_ref/hol_testset.F250_HD_merged.1970"
	log:
		"logs/ref_pop_frqs/F250_HD_merged.txt"
	benchmark:
		"benchmarks/ref_pop_frqs/F250_HD_merged.benchmark.txt"
	output:
		frqx = "ref_vcfs/F250_HD_merged.frqx",
		frq = "ref_vcfs/F250_HD_merged.frq"
	shell: "(plink --bfile {params.iprefix}  --cow --real-ref-alleles --freq -out {params.oprefix};plink --bfile {params.iprefix}  --cow --real-ref-alleles --freqx -out {params.oprefix})>{log}"

rule impref_frq:
	input:
		bed = "merge_ref/hol_testset.F250_HD_merged.1970.bed",
		idslist = "dataprepper/{assay}_ids.list{list}.txt"
	params:
		assayset = "ref_vcfs/{assay}_animals.list{list}",
		iprefix = "merge_ref/hol_testset.F250_HD_merged.1970",
		idslist = " --keep dataprepper/{assay}_ids.list{list}.txt"
	output:
#		assayset = "ref_vcfs/{assay}_animals.list{list}.frq"
	log:
		"logs/downsample_assay_frqs/{assay}_animals.list{list}"
	benchmark:
		"benchmarks/downsample_assay_frqs/{assay}_animals.list{list}.txt"
	shell: "(plink --bfile {params.iprefix}  --cow --real-ref-alleles --freq {params.idslist} -out {params.assayset})>{log}"
	#snakemake -s impacc.snakefile ref_vcfs/f250_animals.list1.frq



rule downsample_assay_frqs:
	input:
		bed = "merge_ref/hol_testset.F250_HD_merged.1970.bed",
		idslist = "dataprepper/{assay}_ids.list{list}.txt"
	params:
		assayset = "ref_vcfs/{assay}_animals.list{list}",
		iprefix = "merge_ref/hol_testset.F250_HD_merged.1970",
		idslist = " --keep dataprepper/{assay}_ids.list{list}.txt"
	output:
		assayset = "ref_vcfs/{assay}_animals.list{list}.frq"
	log:
		"logs/downsample_assay_frqs/{assay}_animals.list{list}"
	benchmark:
		"benchmarks/downsample_assay_frqs/{assay}_animals.list{list}.txt"
	shell: "(plink --bfile {params.iprefix}  --cow --real-ref-alleles --freq {params.idslist} -out {params.assayset})>{log}"
	#snakemake -s impacc.snakefile ref_vcfs/f250_animals.list1.frq


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

rule minimac_decompress:
	input:
		vcf = "minimac_imputed/run{run}/{sample}.chr{chr}.imputed.dose.vcf.gz"
	output:
		vcf = temp("minimac_imputed/run{run}/{sample}.chr{chr}.imputed.dose.vcf")
	shell:
		"gunzip -c {input.vcf} > {output.vcf}"

rule imp_acc:
	input:
		true = "ref_vcfs/F250_HD_merged.chr{chr}.pickle",
		imputed = vcffinder,
		#imputed = "impute2_vcf/run{run}/{sample}.chr{chr}.imputed.vcf",
		frq = "ref_vcfs/F250_HD_merged.chr{chr}.frq",
	params:
		chrom = "{chr}",
		acc = "imp_acc/run{run}/{sample}.txt",
		script = impaccscript
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
		"(python {params.script} {input.true} {input.imputed} {params.acc} {output.corrs})>{log}"


# rule minimac_imp_acc:
# 	input:
# 		true = "ref_vcfs/F250_HD_merged.chr{chr}.pickle",
# 		imputed = "minimac_imputed/run{run}/{sample}.chr{chr}.imputed.dose.vcf",
# 		#imputed = "impute2_vcf/{sample}.chr{chr}.imputed.vcf",
# 		frq = "ref_vcfs/F250_HD_merged.chr{chr}.frq",
# 	params:
# 		chrom = "{chr}",
# 		acc = "imp_acc/run{run}/{sample}.txt",
# 		#vcf = "minimac_imputed/run{run}/{sample}.chr{chr}.imputed.dose.vcf"
# 	log:
# 		"logs/imp_acc/run{run}/{sample}.chr{chr}.txt"
# 	benchmark:
# 		"benchmarks/imp_acc/run{run}/{sample}.chr{chr}.benchmark.txt"
# 	output:
# 		corrs = "imp_acc/run{run}/{sample}.chr{chr}.snp_correlations.csv", # This will contain all of the correlations for each base pair of the assay/run/chromosome
# 	shell:
# 		"(python bin/minimac_impacc.py {input.true} {input.imputed} {params.acc} {output.corrs})>{log}"

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
