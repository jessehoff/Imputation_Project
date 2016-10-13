shell.prefix("set -euo pipefail;")

rule filter_variants_snp50:
	input:
		map="maps/9913_SNP50.map"
	params: 
		genprefix="raw_genotypes/58336.160906.100.test_{sample}",
		outprefix="allele_filtered/{sample}"
	output:
		bed="allele_filtered/{sample}.bed"
	shell:
		"plink --file {params.genprefix} --map {input.map} --cow --not-chr 0 --geno .05 --make-bed --out  {params.outprefix}"

rule filter_variants_ggpf250:
	input:
		map="maps/9913_GGPF250.map"
	params:
		genprefix = "raw_genotypes/227234.160906.100.test_ggpf250_A",
		outprefix = "allele_filtered/ggpf250"
	output:
		bed="allele_filtered/ggpf250.bed"
	shell:
		"plink --file {params.genprefix} --map {input.map} --cow --not-chr 0 --geno .05 --make-bed --out {params.outprefix}"

rule filter_variants_hd:
	input:
                map="maps/9913_HD.map"
	params:
                genprefix = "raw_genotypes/777962.160906.100.test_hd_A",
                outprefix = "allele_filtered/hd"
	output:
                bed="allele_filtered/hd.bed"
	shell:
                "plink --file {params.genprefix} --map {input.map} --cow --not-chr 0 --geno .05 --make-bed --out {params.outprefix}"

rule filter_individuals:
        input:
        	bed="allele_filtered/{sample}.bed"
	params:
                prefix="allele_filtered/{sample}",
                oprefix="individual_filtered/{sample}"
	output:
		bed="individual_filtered/{sample}.bed",
		log="individual_filtered/{sample}.log"
	shell:
                "plink --bfile {params.prefix}  --cow --mind .05 --make-bed --out  {params.oprefix}"
rule report:
	input:
		bed="individual_filtered/{sample}.log"

rule filter_hwe_variants:
	input:
		bed="individual_filtered/{sample}.bed"
	params: 
		prefix="individual_filtered/{sample}",
		oprefix="hwe_filtered/{sample}"
	output:
		bed="hwe_filtered/{sample}.bed"
	shell:
		"plink --bfile {params.prefix} --cow --nonfounders --hwe 0.01 --make-bed --out {params.oprefix}"

