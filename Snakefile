rule filter_variants:
	input:
		map="maps/9913_SNP50.map"
	params: 
		genprefix="raw_genotypes/58336.160906.100.test_snp50_A",

		outprefix="allele_filtered/snp50_a"
	output:
		bed="allele_filtered/snp50_a.bed",
		bim="allele_filtered/snp50_a.bim",
		fam="allele_filtered/snp50_a.fam"
	shell:
		"plink --file {params.genprefix} --map {input.map} --cow --not-chr 0 --geno .05 --make-bed --out  {params.outprefix}"
rule filter_individuals_variants:
        input:
                map="maps/9913_SNP50.map",
        	bed="allele_filtered/snp50_a.bed",
		bim="allele_filtered/snp50_a.bim",
                fam="allele_filtered/snp50_a.fam"
	params:
                prefix="allele_filtered/snp50_a",
                oprefix="individual_filtered/snp50_a"
	output:
                bim="individual_filtered/snp50_a.bim",
                fam="individual_filtered/snp50_a.fam",
		bed="individual_filtered/snp50_a.bed"
	shell:
                "plink --bfile {params.prefix}  --cow --mind .05 --make-bed --out  {params.oprefix}"


