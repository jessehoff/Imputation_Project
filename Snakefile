map_dict = {'777962':"/CIFS/MUG01_N/taylorjerr/PLINK_FILES/9914_HD_161214.map",'227234':"/CIFS/MUG01_N/taylorjerr/PLINK_FILES/9913_GGPF250_161214.map",'58336':"/CIFS/MUG01_N/taylorjerr/PLINK_FILES/9913_SNP50_161214.map", '139977':"/CIFS/MUG01_N/taylorjerr/PLINK_FILES/9913_GGPHDv3_161214.map", '26504':"/CIFS/MUG01_N/taylorjerr/PLINK_FILES/9913_GGPLDv3_161214.map", '30105':"/CIFS/MUG01_N/taylorjerr/PLINK_FILES/9913_GGPLDV4_161214.map",'76999':"/CIFS/MUG01_N/taylorjerr/PLINK_FILES/9913_GGP90KT_161214.map"}
#This map dictionary should be able to remain the same, and we can add new maps for whichever new assays become available in future datasets

#should be adapted to just parse the filename for the info, and have as many dict keys as maps available
def mapdicter(shoein): 
	t = shoein.sample
	num_sites = t.split('.')[0] #the file naming involves the number of sites in the genotype query, and so indicates which map to use. 
	return map_dict[num_sites]

def sampdicter(wildcards):         
	return samp_dict[wildcards.sample]

#This step produces a new ped file that does not have any duplicate individuals.  
#Finds duplicated individuals in ID file.  Only adds these individuals once when parsing ped file line by line.
#Python Script (duplicate_filter.py) takes beginning ID and ped files.  Produces list of duplicate IDs, and then a new ped file both in specified locations
#Input File Locations: raw_genotypes/
#Output File Location: duplicates_filtered/ and duplicates_filtered/dup_ids/
#Output File Types: (.ped, .txt)
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
		"python bin/duplicate_filter.py {input.id} {output.dup} {input.ped} {output.ped}"

#Python Script (duplicate_logging.py) reads the list of duplicated ids produced in filter_duplicate_individuals rule, counts lines, and writes to specified csv that will serve as logging for the rest of the filtering process. 
#Input File Location: raw_genotypes/dup_ids/
#Output File Location: filter_logs/
#output File Type: (.csv)
rule create_filter_log:
        input:
                ped="duplicates_filtered/{sample}.ped",
                log = "duplicates_filtered/dup_ids/{sample}.txt"
        output:
                csv = "filter_logs/{sample}.csv"
        shell:
                "python bin/duplicate_logging.py {input.log} {output.csv}"

#Dictionaries: This rule takes .map and .ped inputs from dictionaries based on specified output 'wildcard' which is the name of the assay.  
#PLINK command (freq) -- Determines call rates for each variant in dataset and returns how many times each allele appears in every assay
#Python Script (allele_call_rate viusalization.py) and (allele_removal_prediction.py).  Takes argument and identifies how many alleles will be removed given a specific filter call rate.  Visualization script creates histograms with distributions of call rates and exports to png file in allele_stats/ directory
#Input File Location: raw_genotypes/
#Output File Location: allele_stats/
#Output File Types: (.frq, .log, .nosex)
rule variant_stats:
	input:
		map = mapdicter,
		dup = "duplicates_filtered/dup_ids/{sample}.txt",
		input = "duplicates_filtered/{sample}.ped",
		log = "filter_logs/{sample}.csv",
		txt = "input_test.txt"
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

#Python Script for Visualization (allele_call_rate_visualization.py) reads in .frq file generated by PLINK's --freq function
#Creates a histogram of the missing genotype rate for each SNP in the dataset.  This call rate is calculated for each SNP as (NCHROBS)/max(NCHROBS)
#Writes figure to a .png file
#Histogram uses 100 bins, and a y-axis max of 5000 SNPs, which look good for existing assays, but can be easily adjusted.
#Input File Location: allele_stats/
#Output File Location: allele_stats/figures/
#Output File Type: (.png)

rule allele_call_rate_visualization:
	input:
		frq = "allele_stats/{sample}.frq"
	output:
		png = "allele_stats/figures/{sample}.individual_call_rate.png"
	shell:
		"python allele_call_rate_visualization.py {input.frq} {output.png}"


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
		stats = "allele_stats/{sample}.frq",
		csv = "filter_logs/{sample}.csv",
		png = "allele_stats/figures/{sample}.individual_call_rate.png"
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
		"plink --file {params.inprefix} --map {input.map} --keep-allele-order --cow --not-chr 0 --exclude maps/map_issues/snp_ids_to_exclude.txt --geno .05 --make-bed --out {params.oprefix}; python bin/allele_filtered_log_parsing.py {output.log} {input.csv}"




#PLINK statistic (missing) -- Gathers statistics on individual call rates (reports proportion of missing genotypes for each individual)
#Input File Location: /allele_filtered/
#Output: output file suffixes will be (.imiss, .lmiss, .log, .nosex)
rule individual_stats:
	input:
		rules.create_filter_log.output,
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
		"plink --bfile {params.inprefix} --cow --missing --keep-allele-order --out {params.oprefix}" #python bin/individual_call_rate_visualization.py"

#Python Script for Visualization (individual_call_rate_visualization.py) reads in .imiss file generated by PLINK's --missing function
#Creates a histogram of the missing genotyping rate for each individual in the dataset (F_MISS)
#Writes figure to a .png file
#Histogram uses 50 bins, and a y-axis max of 50 animals, which look good for existing assays, but can be easily adjusted.
#Input File Location: individual_stats/
#Output File Location: individual_stats/figures/
#Output File Type: (.png)
rule individual_call_rate_visualization:
	input:
		imiss="individual_stats/{sample}.imiss"
	output:
		png="individual_stats/figures/{sample}.individual_call_rate.png"
	shell:
		"python bin/individual_call_rate_visualization.py {input.imiss} {output.png}"

#PLINK filter (mind) -- Filters individuals based on specified call rate.  Individuals missing more than given proportion of their genotypes will be filtered out of the dataset.
#Python script (individual_filtered_log_parsing.py) -- parses output logfile for important metadata and saves to csv in /filter_logs directory.
#Input File Location: /allele_filtered/
#Output File Location: /individual_filtered
#Output: output file suffixes will be (.bed, .bim, .fam, .irem, .log, .nosex)

rule filter_individuals:
        input:
        	rules.create_filter_log.output,
		csv="filter_logs/{sample}.csv",
		bed="allele_filtered/{sample}.bed",
		imiss="individual_stats/{sample}.imiss",
		png="individual_stats/figures/{sample}.individual_call_rate.png"
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
		"plink --bfile {params.inprefix} --cow --mind .05 --keep-allele-order --make-bed --out {params.oprefix}; python bin/individual_filtered_log_parsing.py {output.log} {input.csv}"

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

rule hwe_visualization:
	input:
		hwe="hwe_stats/{sample}.hwe"	
	output:
		png="hwe_stats/figures/{sample}.hwe_pvalues.png"
	shell:
		"python bin/hwe_visualization.py {input.hwe} {output.png}"
#PLINK filter (hwe) -- Filters loci based on their HWE P values. Variants with HWE P values below a specified value will be removed from the dataset
#Python script (hwe_log_parsing.py) -- parses output logfile for important metadata and saves to csv in /filter_logs directory
#Input File Location: /individual_filtered/
#Output File Location: /hwe_filtered/
#Output: output file suffixes will be (.bed, .bim, .fam, .irem, .log, .nosex)
rule filter_hwe_variants:
	input:
		bed="individual_filtered/{sample}.bed",
		stats="hwe_stats/{sample}.hwe",
		csv="filter_logs/{sample}.csv",
		png="hwe_stats/figures/{sample}.hwe_pvalues.png"
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
		"plink --bfile {params.inprefix} --cow --nonfounders  --keep-allele-order --hwe 0.00000001 --make-bed --out {params.oprefix}; python bin/hwe_log_parsing.py {output.log} {input.csv}"

#PLINK function (impute-sex) looks at sex provided by ped file, and at X chromosome heterozygosity (and y chromosome variants if provided), and determines whether an animal is male, female, or unknown. If sex from ped file is unknown, this will impute the sex if possible, and write that into the new bed file that it produces. 
#Python Scripts: (missexed_animals_filter.py) reads the .sexcheck file produced by impute sex function.  If an individual's sex is changed from known M/F to the opposite sex, it's ID will be written to the {sample}.missexed_animals.txt file, and removed in subsequent step.  (missexed_animals_filter_logging.py) will count lines of filter output, and write to the assay's csv log file.  
##Input File Location: hwe_filtered/
##Output File Location: sex_impute/
##Output: output file suffixes will be (.bed, .bim, .fam, .log, .sexcheck, .missexed_animals.txt)
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
		"plink --bfile {params.inprefix} --cow --impute-sex ycount --make-bed --out {params.oprefix}; python bin/missexed_animals_filter.py {output.sexcheck} {output.txt}; python bin/missexed_animals_filter_logging.py {output.txt} {input.csv}" #{params.logprefix}.csv"

#PLINK function (remove) takes "sex_impute/{sample}.missexed_animals.txt" and produces new bed file without listed IDs, removing missexed animals
#Input File Locations: sex_impute/
#Output File Locations: correct_sex/
#Output File Types: (.bed, .bim, .fam, .log, .hh, .nosex)

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


#This list contains the different file names that will be merged together.  
DATA =['139977.170112.2326.GGPHDV3', '227234.170112.325.GGPF250', '26504.170112.3126.GGPLDV3', '30105.170112.2500.GGPLDV4', '58336.170112.315.SNP50A', '58336.170112.335.SNP50B', '58336.170112.3399.SNP50C', '76999.170112.3498.GGP90KT']

#Python Script (file_list_maker.py) searches the correct_sex sub-directory for all of the files made with a .bed suffix, and then writes the file names to a .txt file (allfiles.txt)
#PLINK function (merge-list) takes this list as an input and merges the different files into a single .bed file
#Input File Locations: correct_sex/
#Output File Location: merged_files/
#Output File Types: (.bed, .bim, .fam, .log, .hh, .nosex)
rule merge_assays:
	input:
		expand("correct_sex/{previous}.bed", previous=DATA), #Expand function allows us to ask for one file -- ___merged.bed, and have it produce all of the necessary prerequesite files from the DATA list above.
		expand("correct_sex/{previous}.bim", previous=DATA),
		expand("correct_sex/{previous}.fam", previous=DATA),
		expand("correct_sex/{previous}.log", previous=DATA)
	params:
		oprefix="merged_files/170112_merged"	
	output:
		#mergefilelist= "correct_sex/allfiles.txt"
                bim = "merged_files/170112_merged.bim",
                fam = "merged_files/170112_merged.fam",
                log = "merged_files/170112_merged.log",
                bed = "merged_files/170112_merged.bed"
	shell:
		"python file_list_maker.py correct_sex/allfiles.txt; plink --merge-list correct_sex/allfiles.txt  --cow --make-bed --out {params.oprefix}"



#Chrsplit -- Splits final .bed output of PLINK filtering septs into individual chromosomes to prepare for phasing step
#Input file location: merged_files/
#Output file location: chrsplit/
#Output file types: (.bed, .bim, .fam, .log, .hh, .nosex)
rule split_chromosomes:
	input:
		bed = "merged_files/170112_merged.bed",
		bim = "merged_files/170112_merged.bim",
		fam = "merged_files/170112_merged.fam",
		log = "merged_files/170112_merged.log"
	params:
		inprefix = "merged_files/170112_merged",
		oprefix = "chrsplit/170112_merged.chr"
	benchmark:                 
		"filter_benchmarks/split_chromosomes170112_merged.txt"
	output:
		bed = "chrsplit/170112_merged.chr1.bed", #may need to make this an oprefix param
		bim = "chrsplit/170112_merged.chr1.bim",
		fam = "chrsplit/170112_merged.chr1.fam",
		log = "chrsplit/170112_merged.chr1.log"
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
		

















		 
