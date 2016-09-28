rule filter_variants:
	input:
		map="maps/9913_SNP50.map"
	params: 
		genprefix="raw_genotypes/58336.160906.100.test_{sample}",

		outprefix="allele_filtered/{sample}"
	output:
		bed="allele_filtered/{sample}.bed",
		bim="allele_filtered/{sample}.bim",
		fam="allele_filtered/{sample}.fam"
	shell:
		"plink --file {params.genprefix} --map {input.map} --cow --not-chr 0 --geno .05 --make-bed --out  {params.outprefix}"
rule filter_individuals_variants:
        input:
                map="maps/9913_SNP50.map",
        	bed="allele_filtered/{sample}.bed",
		bim="allele_filtered/{sample}.bim",
                fam="allele_filtered/{sample}.fam"
	params:
                prefix="allele_filtered/{sample}",
                oprefix="individual_filtered/{sample}"
	output:
                bim="individual_filtered/{sample}.bim",
                fam="individual_filtered/{sample}.fam",
		bed="individual_filtered/{sample}.bed"
	shell:
                "plink --bfile {params.prefix}  --cow --mind .05 --make-bed --out  {params.oprefix}"


