map_dict = {'777962':"/CIFS/MUG01_N/taylorjerr/PLINK_FILES/9914_HD_161214.map",'227234':"/CIFS/MUG01_N/taylorjerr/PLINK_FILES/9913_GGPF250_161214.map",'58336':"/CIFS/MUG01_N/taylorjerr/PLINK_FILES/9913_SNP50_161214.map", '139977':"/CIFS/MUG01_N/taylorjerr/PLINK_FILES/9913_GGPHDv3_161214.map", '26504':"/CIFS/MUG01_N/taylorjerr/PLINK_FILES/9913_GGPLDv3_161214.map", '30105':"/CIFS/MUG01_N/taylorjerr/PLINK_FILES/9913_GGPLDV4_161214.map",'76999':"/CIFS/MUG01_N/taylorjerr/PLINK_FILES/9913_GGP90KT_161214.map"}


#should be adapted to just parse the filename for the info, and have as many dict keys as maps available
def mapdicter(shoein): 
	t = shoein.sample
	num_sites = t.split('.')[0] #the file naming involves the number of sites in the genotype query, and so indicates which map to use. 
	return map_dict[num_sites]

def sampdicter(wildcards):         
	return samp_dict[wildcards.sample]


rule filter_duplicate_individuals:
	input:
		id = "raw_genotypes/{sample}.ID",
		ped = "raw_genotypes/{sample}.ped"
	benchmark:
		"filter_benchmarks/duplicates_filtered/{sample}.txt"
	output:
		dup = "duplicates_filtered/dup_ids/{sample}.txt",
		ped = "duplicates_filtered/{sample}.ped"
	shell:
		"python duplicates_filtered/duplicate_filter.py {input.id} {output.dup} {input.ped} {output.ped}"

#Dictionaries: This rule takes .map and .ped inputs from dictionaries based on specified output 'wildcard' which is the name of the assay.  
#PLINK command (freq) -- Determines call rates for each variant in dataset and returns how many times each allele appears in every assay
#Python Script (allele_call_rate viusalization.py) and (allele_removal_prediction.py).  Takes argument and identifies how many alleles will be removed given a specific filter call rate.  Visualization script creates histograms with distributions of call rates and exports to png file in allele_stats/ directory
#Input File Location: raw_genotypes/
#Output File Location: allele_stats/
#Output File Types: (.frq, .log, .nosex)

#rule mergeraws: #adds a mix of sample files. Ma
#	--keep ./checkflip/snp50_a.fam


rule variant_stats:
	input:
		map = mapdicter,
		dup = "duplicates_filtered/dup_ids/{sample}.txt",
		input = "duplicates_filtered/{sample}.ped"
	params:
		inprefix = "duplicates_filtered/{sample}",
		oprefix = "allele_stats/{sample}"
	benchmark:
		"filter_benchmarks/variant_stats/{sample}.txt"
	output: 
		frq = "allele_stats/{sample}.frq",
		log = "allele_stats/{sample}.log",
		hh = "allele_stats/{sample}.hh"
	shell:
		"plink --file {params.inprefix} --map {input.map} --keep-allele-order --cow --nonfounders --freq --out {params.oprefix}" 


#rule allele_call_rate_visualization:
#	input:
#		frq = "allele_stats/{sample}.frq"
#	output:
#		png = "allele_stats/{sample}.png"
#	shell:
#		"python allele_call_rate_visualization.py {input.frq}"



#Dictionaries for selecting correct .map and .ped files:  This uses two functions, mapdicter and sampdicter to match the correct .ped and .map files with each assay. Using a specific 'wildcard' for the output of this rule will call the correct input map and ped file for the filter to be used on. 
#PLINK filter (mind) -- Filters individuals based on specified call rate.  Individuals missing more than given proportion of their genotypes will be filtered out of the dataset. 
#Python script (allele_filtered_log_parsing.py) -- parses output logfile for important metadata and saves to assay's specific csv in /filter_logs directory. 
#Input File Location: /raw_genotypes/ 
#Output File Location: /allele_filtered 
#Output: output file suffixes will be (.bed, .bim, .fam, .irem, .log, .nosex)

rule filter_variants:
	input:	
		map = mapdicter,
		dup = "duplicates_filtered/dup_ids/{sample}.txt",
                input = "duplicates_filtered/{sample}.ped",
		stats = "allele_stats/{sample}.frq"
	params:
		inprefix = "duplicates_filtered/{sample}",
		oprefix="allele_filtered/{sample}",
		logprefix="filter_logs/{sample}"
	benchmark:                 
		"filter_benchmarks/filter_variants/{sample}.txt"
	output:
		bed="allele_filtered/{sample}.bed",
		bim="allele_filtered/{sample}.bim",
		fam="allele_filtered/{sample}.fam",
		log="allele_filtered/{sample}.log"
	shell:
		"plink --file {params.inprefix} --map {input.map} --keep-allele-order --cow --not-chr 0 --exclude maps/map_issues/snp_ids_to_exclude.txt --geno .05 --make-bed --out  {params.oprefix}"



rule variant_filter_log:
	input:
		bed="allele_filtered/{sample}.bed",
		log = "allele_filtered/{sample}.log"
	output:
		csv = "filter_logs/{sample}.csv"
	shell:
		"python allele_filtered/allele_filtered_log_parsing.py {input.log} {output.csv}"



#PLINK statistic (missing) -- Gathers statistics on individual call rates (reports proportion of missing genotypes for each individual)
#Input File Location: /allele_filtered/
#Output: output file suffixes will be (.imiss, .lmiss, .log, .nosex)
rule individual_stats:
	input:
		rules.variant_filter_log.output,
		bed="allele_filtered/{sample}.bed"
	params:
		inprefix="allele_filtered/{sample}",
		oprefix="individual_stats/{sample}"
	benchmark:                 
		"filter_benchmarks/individual_stats/{sample}.txt"
	output:
		imiss="individual_stats/{sample}.imiss",
		lmiss="individual_stats/{sample}.lmiss",
		log ="individual_stats/{sample}.log"
	shell:
		"plink --bfile {params.inprefix} --cow --missing --keep-allele-order --out {params.oprefix}" #python individual_call_rate_visualization.py"

rule individual_call_rate_visualization:
	input:
		imiss="individual_stats/{sample}.imiss"
	output:
		png="individual_stats/{sample}.png"
	shell:
		"python individual_call_rate_visualization.py"

#PLINK filter (mind) -- Filters individuals based on specified call rate.  Individuals missing more than given proportion of their genotypes will be filtered out of the dataset.
#Python script (individual_filtered_log_parsing.py) -- parses output logfile for important metadata and saves to csv in /filter_logs directory.
#Input File Location: /allele_filtered/
#Output File Location: /individual_filtered
#Output: output file suffixes will be (.bed, .bim, .fam, .irem, .log, .nosex)

rule filter_individuals:
        input:
        	rules.variant_filter_log.output,
		csv="filter_logs/{sample}.csv",
		bed="allele_filtered/{sample}.bed",
		imiss="individual_stats/{sample}.imiss"
	params:
                inprefix="allele_filtered/{sample}",
                oprefix="individual_filtered/{sample}",
	benchmark:                 
		"filter_benchmarks/filter_individuals/{sample}.txt"
	output:
		bed="individual_filtered/{sample}.bed",
		bim="individual_filtered/{sample}.bim",
		fam="individual_filtered/{sample}.fam",
		log="individual_filtered/{sample}.log"
		
	shell:
		"plink --bfile {params.inprefix} --cow --mind .05 --keep-allele-order --make-bed --out {params.oprefix}; python individual_filtered/individual_filtered_log_parsing.py {output.log} {input.csv}"



#Gathers statistics on individual Hardy Weinberg Equilibrium (reports HWE P value at each locus)
#Input File Location: /individual_filtered/
#Output File Location: /hwe_stats/
#Output: output file suffixes will be (.hwe, .log, .nosex)

rule hwe_stats:
	input:
		csv="filter_logs/{sample}.csv",
		bed="individual_filtered/{sample}.bed"
	params:
		inprefix="individual_filtered/{sample}",
		oprefix="hwe_stats/{sample}"
	benchmark:                 
		"filter_benchmarks/hwe_stats/{sample}.txt"
	output:
		hwe="hwe_stats/{sample}.hwe",
		log="hwe_stats/{sample}.log"
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
		stats="hwe_stats/{sample}.hwe",
		csv="filter_logs/{sample}.csv"
	params: 
		inprefix="individual_filtered/{sample}",
		oprefix="hwe_filtered/{sample}"
	benchmark:                 
		"filter_benchmarks/filter_hwe_variants/{sample}.txt"
	output:
		bed="hwe_filtered/{sample}.bed",
		bim="hwe_filtered/{sample}.bim",
		fam="hwe_filtered/{sample}.fam",
		log="hwe_filtered/{sample}.log"
	shell:
		"plink --bfile {params.inprefix} --cow --nonfounders  --keep-allele-order --hwe 0.00000001 --make-bed --out {params.oprefix}; python hwe_filtered/hwe_log_parsing.py {output.log} {input.csv}"


rule impute_sex:
	input:
		bed="hwe_filtered/{sample}.bed",
		csv="filter_logs/{sample}.csv"
	params:
		logprefix="filter_logs/{sample}",
		inprefix="hwe_filtered/{sample}",
		oprefix="sex_impute/{sample}"
	benchmark:
		"filter_benchmarks/impute_sex/{sample}.txt"
	output:
		bed="sex_impute/{sample}.bed",
		bim="sex_impute/{sample}.bim",
		fam="sex_impute/{sample}.fam",
		log="sex_impute/{sample}.log",
		sexcheck="sex_impute/{sample}.sexcheck",
		txt="sex_impute/{sample}.missexed_animals.txt"
	shell:
		"plink --bfile {params.inprefix} --cow --impute-sex ycount --make-bed --out {params.oprefix}; python sex_impute/missexed_animals_filter.py {output.sexcheck} {output.txt}; python sex_impute/missexed_animals_filter_logging.py {output.txt} {input.csv}" #{params.logprefix}.csv"



rule remove_missexed_animals:
	input:
		txt="sex_impute/{sample}.missexed_animals.txt",
		bed="sex_impute/{sample}.bed"
	params:
		inprefix="sex_impute/{sample}",
		oprefix="correct_sex/{sample}"
	output:
		bed="correct_sex/{sample}.bed",
                bim="correct_sex/{sample}.bim",
                fam="correct_sex/{sample}.fam",
                log="correct_sex/{sample}.log",	
	shell:
		"plink --bfile {params.inprefix} --cow --remove {input.txt} --make-bed --out {params.oprefix}"


#Mendel Error Rates will happen last before merging across assays

DATA = ['58336.161214.315.SNP50A', '58336.161214.335.SNP50B', '58336.161214.3399.SNP50C'] 
#DATA =['139977.161214.2326.GGPHDV3', '227234.161214.325.GGPF250', '26504.161214.3126.GGPLDV3', '30105.161214.2500.GGPLDV4', '58336.161214.315.SNP50A', '58336.161214.335.SNP50B', '58336.161214.3399.SNP50C', '76999.161214.3498.GGP90KT']
#DATA = ['58336.161214.335.SNP50B', '227234.161214.325.GGPF250']
rule merge_assays: #split this into two steps
	input:
		expand("correct_sex/{previous}.bed", previous=DATA),
		expand("correct_sex/{previous}.bim", previous=DATA),
		expand("correct_sex/{previous}.fam", previous=DATA),
		expand("correct_sex/{previous}.log", previous=DATA)
	params:
		oprefix="merged_files/161214_merged"	
#	oprefix="merged_files/{merged}"
#	benchmark:                 
#		"filter_benchmarks/merge_assays/{merged}.txt"
	output:
		#mergefilelist= "correct_sex/allfiles.txt"
                bim = "merged_files/161214_merged.bim",
                fam = "merged_files/161214_merged.fam",
                log = "merged_files/161214_merged.log",
                bed = "merged_files/161214_merged.bed"
	shell:
		"python file_list_maker.py correct_sex/allfiles.txt; plink --merge-list correct_sex/allfiles.txt  --cow --make-bed --out {params.oprefix}"



#Chrsplit -- Splits final .bed output of PLINK filtering septs into individual chromosomes
#Input file location: ./hwe_filtered/
#Input file types: (.bed, .bim, .fam)
#Output file location: ./shapetest/
#Output file types: (.phased.haps, .phased.sample)
rule split_chromosomes:
	input:
		bed = "merged_files/161214_merged.bed",
		bim = "merged_files/161214_merged.bim",
		fam = "merged_files/161214_merged.fam",
		log = "merged_files/161214_merged.log"
	params:
		inprefix = "merged_files/161214_merged",
		oprefix = "chrsplit/161214_merged.chr"
	benchmark:                 
		"filter_benchmarks/split_chromosomes161214_merged.txt"
	output:
		bed = "chrsplit/161214_merged.chr1.bed", #may need to make this an oprefix param
		bim = "chrsplit/161214_merged.chr1.bim",
		fam = "chrsplit/161214_merged.chr1.fam",
		log = "chrsplit/161214_merged.chr1.log"
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
		

















		 
