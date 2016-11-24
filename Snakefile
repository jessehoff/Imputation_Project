map_dict = {'snp50_a':"maps/9913_SNP50.map", 'snp50_b':"maps/9913_SNP50.map", 'snp50_c':"maps/9913_SNP50.map", 'hd':"maps/9913_HD.map", 'ggpf250':"maps/9913_GGPF250.map" }

samp_dict = {'snp50_a': "raw_genotypes/58336.160906.100.test_snp50_A", 'snp50_b':"raw_genotypes/58336.160906.100.test_snp50_B", 'snp50_c':"raw_genotypes/58336.160906.100.test_snp50_C", 'ggpf250':"raw_genotypes/227234.160906.100.test_ggpf250_A", 'hd':"raw_genotypes/777962.160906.100.test_hd_A"}

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
	benchmark:
		"filter_benchmarks/variant_stats/{sample}.txt"
	output: 
		frq = "allele_stats/{sample}.frq"
	shell:
		"plink --file {params.inprefix} --map {input.map} --keep-allele-order --cow --nonfounders --freq --out {params.oprefix}" # python allele_call_rate_visualization.py"



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
	benchmark:                 
		"filter_benchmarks/filter_variants/{sample}.txt"
	output:
		bed="allele_filtered/{sample}.bed"
	shell:
		"plink --file {params.inprefix} --map {input.map} --keep-allele-order --cow --not-chr 0 --exclude maps/map_issues/snp_ids_to_exclude.txt --geno .05 --make-bed --out  {params.oprefix}; python allele_filtered/allele_filtered_log_parsing.py {params.oprefix}.log {params.logprefix}.csv"



#PLINK statistic (missing) -- Gathers statistics on individual call rates (reports proportion of missing genotypes for each individual)
#Input File Location: /allele_filtered/
#Output: output file suffixes will be (.imiss, .lmiss, .log, .nosex)
rule individual_stats:
	input:
		bed="allele_filtered/{sample}.bed"
	params:
		inprefix="allele_filtered/{sample}",
		oprefix="individual_stats/{sample}"
	benchmark:                 
		"filter_benchmarks/individual_stats/{sample}.txt"
	output:
		imiss="individual_stats/{sample}.imiss"
	shell:
		"plink --bfile {params.inprefix} --cow --missing --keep-allele-order --out {params.oprefix}" #python individual_call_rate_visualization.py"



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
	benchmark:                 
		"filter_benchmarks/filter_individuals/{sample}.txt"
	output:
		bed="individual_filtered/{sample}.bed"
		
	shell:
		"plink --bfile {params.inprefix} --cow --mind .05 --keep-allele-order --make-bed --out {params.oprefix}; python individual_filtered/individual_filtered_log_parsing.py {params.oprefix}.log {params.logprefix}.csv"




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
	benchmark:                 
		"filter_benchmarks/hwe_stats/{sample}.txt"
	output:
		hwe="hwe_stats/{sample}.hwe"
	shell:
		"plink --bfile {params.inprefix}  --keep-allele-order --cow --nonfounders --hardy --out {params.oprefix}"



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
	benchmark:                 
		"filter_benchmarks/filter_hwe_variants/{sample}.txt"
	output:
		bed="hwe_filtered/{sample}.bed"

	shell:
		"plink --bfile {params.inprefix} --cow --nonfounders  --keep-allele-order --hwe 0.01 --make-bed --out {params.oprefix}" #python hwe_filtered/hwe_log_parsing.py {params.oprefix}.log {params.logprefix}.csv"


rule missexed_filter:
	




#Mendel Error Rates will happen last before merging across assays

DATA = ['snp50_a', 'snp50_b', 'snp50_c', 'hd', 'ggpf250']

rule merge_assays:
	input:
		expand("hwe_filtered/{previous}.bed", previous=DATA)
	params:
		oprefix="merged_files/merged"
	benchmark:                 
		"filter_benchmarks/merge_assays/merged.txt"
	output:
		"merged_files/merged.bed"
	shell:
		"python hwe_filtered/file_list_maker.py; plink --merge-list hwe_filtered/allfiles.txt  --cow --make-bed --out {params.oprefix}"







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
	benchmark:                 "filter_benchmarks/split_chromosomes/{sample}.txt"
	output:
		bed = "chrsplit/{sample}.chr1.bed" #may need to make this an oprefix param
	shell:
		"for chr in $(seq 1 30); do plink --bfile {params.inprefix}  --keep-allele-order --chr $chr --make-bed  --nonfounders --cow --out {params.oprefix}$chr; done"









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
		

















		 
