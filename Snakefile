#def Map_Matcher(startingfilename): there will be a single wildcard object, which should correspond to the sample name. The sample name should start with a number so the funciton should use that number to return a map name
#	Mapfiledict = {227234: ['map="maps/9913_ggpf250.map"'],etc etc}
#	return Mapfiledict[startingfilename] this should be the string of the file map, and so the input should 

#def map_matcher(wildcards):
 #  if~ 
	  

#rule variant_statistics:
#	input:


rule filter_variants_snp50:
	input:
		map="maps/9913_SNP50.map"
		#add the freq stats as an input and all other stats as inputs
	params: 
		inprefix="raw_genotypes/58336.160906.100.test_{sample}",
		oprefix="allele_filtered/{sample}"
	output:
		bed="allele_filtered/{sample}.bed"
	shell:
		"plink --file {params.inprefix} --map {input.map} --cow --not-chr 0 --geno .05 --make-bed --out  {params.oprefix}"

rule filter_variants_ggpf250:
	input:
		map="maps/9913_GGPF250.map"
	params:
		inprefix = "raw_genotypes/227234.160906.100.test_ggpf250_A",
		oprefix = "allele_filtered/ggpf250"
	output:
		bed="allele_filtered/ggpf250.bed"
	shell:
		"plink --file {params.inprefix} --map {input.map} --cow --not-chr 0 --geno .05 --make-bed --out {params.oprefix}"

rule filter_variants_hd:
	input:
                map="maps/9913_HD.map" #you will need a function at this point to provide the appropriate map for each samplei
		#Map_matcher
	params:
                inprefix = "raw_genotypes/777962.160906.100.test_hd_A",
                oprefix = "allele_filtered/hd"
	output:
                bed="allele_filtered/hd.bed",
	#this log and all other log output names need to be specific to the sample
	shell:
                "plink --file {params.inprefix} --map {input.map} --cow --not-chr 0 --geno .05 --make-bed --out {params.oprefix}"




#PLINK statistic (missing) -- Gathers statistics on individual call rates (reports proportion of missing genotypes for each individual)
#Input File Location: /allele_filtered/
#Output: output file suffixes will be (.imiss, .lmiss, .log, .nosex)
rule individual_stats:
	input:
		bed="allele_filtered/{sample}.bed"
	params:
		inprefix="allele_filtered/{sample}",
		oprefix="individual_stats/{sample}"
	output:
		imiss="individual_stats/{sample}.imiss"
	shell:
		"plink --bfile {params.inprefix} --cow --missing --out {params.oprefix}"



#PLINK filter (mind) -- Filters individuals based on specified call rate.  Individuals missing more than given proportion of their genotypes will be filtered out of the dataset.
#Python script (individual_filtered_log_parsing.py) -- parses output logfile for important metadata and saves to csv in /filter_logs directory.
#Input File Location: /allele_filtered/
#Output File Location: /individual_filtered
#Output: output file suffixes will be (.bed, .bim, .fam, .irem, .log, .nosex)

rule filter_individuals:
        input:
        	bed="allele_filtered/{sample}.bed",
		imiss="individual_stats/{sample}.imiss"
	params:
                inprefix="allele_filtered/{sample}",
                oprefix="individual_filtered/{sample}",
		logprefix="filter_logs/{sample}"
	output:
		bed="individual_filtered/{sample}.bed"
		
	shell:
		"plink --bfile {params.inprefix} --cow --mind .05 --make-bed --out {params.oprefix}; python individual_filtered/individual_filtered_log_parsing.py {params.oprefix}.log {params.logprefix}.csv"




#Gathers statistics on individual Hardy Weinberg Equilibrium (reports HWE P value at each locus)
#Input File Location: /individual_filtered/
#Output File Location: /hwe_stats/
#Output: output file suffixes will be (.hwe, .log, .nosex)

rule hwe_stats:
	input:
		bed="individual_filtered/{sample}.bed"
	params:
		inprefix="individual_filtered/{sample}",
		oprefix="hwe_stats/{sample}"
	output:
		hwe="hwe_stats/{sample}.hwe"
	shell:
		"plink --bfile {params.inprefix} --cow --nonfounders --hardy --out {params.oprefix}"



#PLINK filter (hwe) -- Filters loci based on their HWE P values. Variants with HWE P values below a specified value will be removed from the dataset
#Python script (hwe_log_parsing.py) -- parses output logfile for important metadata and saves to csv in /filter_logs directory
#Input File Location: /individual_filtered/
#Output File Location: /hwe_filtered/
#Output: output file suffixes will be (.bed, .bim, .fam, .irem, .log, .nosex)
rule filter_hwe_variants:
	input:
		bed="individual_filtered/{sample}.bed"
	params: 
		inprefix="individual_filtered/{sample}",
		oprefix="hwe_filtered/{sample}",
		logprefix="filter_logs/{sample}"
	output:
		bed="hwe_filtered/{sample}.bed"

	shell:
		"plink --bfile {params.inprefix} --cow --nonfounders --hwe 0.01 --make-bed --out {params.oprefix}; python hwe_filtered/hwe_log_parsing.py {params.oprefix}.log {params.logprefix}.csv"

