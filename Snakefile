map_dict = {'snp50_a':"maps/9913_SNP50.map", 'snp50_b':"maps/9913_SNP50.map", 'snp50_c':"maps/9913_SNP50.map", 'hd':"maps/9913_HD.map", 'ggpf250':"maps/9913_GGPF250.map" }

samp_dict = {'snp50_a': "raw_genotypes/58336.160906.100.test_snp50_A", 'snp50_b':"raw_genotypes/58336.160906.100.test_snp50_B", 'snp50_c':"raw_genotypes/58336.160906.100.test_snp50_C", 'ggpf250':"raw_genotypes/227234.160906.100.test_ggpf250_A", 'hd':"777962.160906.100.test_hd_A"}

def mapdicter(shoein): 
	return map_dict[shoein.sample]

def sampdicter(wildcards):         
	return samp_dict[wildcards.sample]


#Dictionaries: This rule takes .map and .ped inputs from dictionaries based on specified output 'wildcard' which is the name of the assay.  
#PLINK command (freq) -- Determines call rates for each variant in dataset and returns how many times each allele appears in every assay
#Python Script (allele_call_rate viusalization.py) and (allele_removal_prediction.py).  Takes argument and identifies how many alleles will be removed given a specific filter call rate.  Visualization script creates histograms with distributions of call rates and exports to png file in allele_stats/ directory
#Input File Location: raw_genotypes/
#Output File Location: allele_stats/
#Output File Types: (.frq, .log, .nosex)

rule variant_stats:
	input:
		map = mapdicter
	params:
		inprefix = sampdicter,
		oprefix = "allele_stats/{sample}"
	output: 
		frq = "allele_stats/{sample}.frq"
	shell:
		"plink --file {params.inprefix} --map {input.map} --cow --nonfounders --freq --out {params.oprefix}" # python allele_call_rate_visualization.py"



#Dictionaries for selecting correct .map and .ped files:  This uses two functions, mapdicter and sampdicter to match the correct .ped and .map files with each assay. Using a specific 'wildcard' for the output of this rule will call the correct input map and ped file for the filter to be used on. 
#PLINK filter (mind) -- Filters individuals based on specified call rate.  Individuals missing more than given proportion of their genotypes will be filtered out of the dataset. 
#Python script (allele_filtered_log_parsing.py) -- parses output logfile for important metadata and saves to assay's specific csv in /filter_logs directory. 
#Input File Location: /raw_genotypes/ 
#Output File Location: /allele_filtered 
#Output: output file suffixes will be (.bed, .bim, .fam, .irem, .log, .nosex)
rule filter_variants:
	input:	
		map = mapdicter,
		stats = "allele_stats/{sample}.frq"
	params:
		inprefix = sampdicter,
		oprefix="allele_filtered/{sample}",
		logprefix="filter_logs/{sample}"
	output:
		bed="allele_filtered/{sample}.bed"
	shell:
		"plink --file {params.inprefix} --map {input.map} --cow --not-chr 0 --geno .05 --make-bed --out  {params.oprefix}; python allele_filtered/allele_filtered_log_parsing.py {params.oprefix}.log {params.logprefix}.csv"



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
		"plink --bfile {params.inprefix} --cow --missing --out {params.oprefix}" #python individual_call_rate_visualization.py"



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
		bed="individual_filtered/{sample}.bed",
		stats="hwe_stats/{sample}.hwe"
	params: 
		inprefix="individual_filtered/{sample}",
		oprefix="hwe_filtered/{sample}",
		logprefix="filter_logs/{sample}"
	output:
		bed="hwe_filtered/{sample}.bed"

	shell:
		"plink --bfile {params.inprefix} --cow --nonfounders --hwe 0.01 --make-bed --out {params.oprefix}; python hwe_filtered/hwe_log_parsing.py {params.oprefix}.log {params.logprefix}.csv"


rule missexed_filter:
	




#Mendel Error Rates will happen last before merging across assays



rule merge_assays:
	input:
		bed = "hwe_filtered/{sample}.bed",
		bim = "hwe_filtered/{sample}.bim",
		fam = "hwe_filtered/{sample}.fam"
	output:
		bed = "merged_files/merged.bed"
	shell:
		"plink --bfile assay[0] --bmerge  --cow --make-bed --out merged"







#Chrsplit -- Splits final .bed output of PLINK filtering septs into individual chromosomes
#Input file location: ./hwe_filtered/
#Input file types: (.bed, .bim, .fam)
#Output file location: ./shapetest/
#Output file types: (.phased.haps, .phased.sample)
rule split_chromosomes:
	input:
		bed = "hwe_filtered/{sample}.bed"
	params:
		inprefix = "hwe_filtered/{sample}",
		oprefix = "chrsplit/{sample}.chr"
	output:
		bed = "chrsplit/{sample}.chr1.bed" #may need to make this an oprefix param
	shell:
		"for chr in $(seq 1 30); do plink --bfile {params.inprefix} --chr $chr --make-bed  --nonfounders --cow --out {params.oprefix}$chr; done"









rule run_shapeit:
	input:
		bed = "chrsplit/{sample}.chr1.bed",
		#bim = "chrsplit/{sample}.chr1.bim",
		#fam = "chrsplit/{sample}.chr1.fam"
	params:
		inprefix = "chrsplit/{sample}.chr",
		oprefix="shapetest/{sample}.chr"
	output:
		haps = "shapetest/{sample}.chr1.phased.haps",
		sample = "shapetest/{sample}.chr1.phased.sample"
	shell:
		"for chr in $(seq 1 30); do shapeit --input-bed {params.inprefix}$chr.bed {params.inprefix}$chr.bim {params.inprefix}$chr.fam --duohmm --output-max {params.oprefix}$chr.phased.haps {params.oprefix}$chr.phased.sample; done"
		

















		 
