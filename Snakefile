SNP50 = ['58336.170626.1009.104.A', '58336.170626.108.124.A', '58336.170626.1082.100.C', '58336.170626.11.112.A', '58336.170626.1113.112.C', '58336.170626.1241.101.C', '58336.170626.127.100.B', '58336.170626.1308.103.C', '58336.170626.13134.200.Z', '58336.170626.1339.200.A', '58336.170626.15.129.B', '58336.170626.16709.200.X', '58336.170626.191.102.C', '58336.170626.2090.550.C', '58336.170626.2203.101.A', '58336.170626.273.550.B', '58336.170626.3.104.B', '58336.170626.3.112.B', '58336.170626.32.101.B', '58336.170626.3245.104.C', '58336.170626.395.129.C', '58336.170626.42.129.A', '58336.170626.4359.200.Y', '58336.170626.45.102.B', '58336.170626.67.124.C', '58336.170626.770.103.A', '58336.170626.7744.100.A', '58336.170626.85.550.A']
#'58336.170626.2714.102.A' having issues with ped file in filtering
#'58336.170626.46.103.B', having issues with no sex chromosomes
REF = ['227234.170619.1255.129_A', '227234.170619.12703.100_A', '227234.170619.1351.101_A', '227234.170619.1667.550_A', '227234.170619.172.103_A', '227234.170619.1994.200_A', '227234.170619.219.102_A', '227234.170619.442.112_A', '227234.170619.500.104_A', '227234.170619.74.124_A', '777962.170619.10.129_C', '777962.170619.11.550_A', '777962.170619.136.124_A', '777962.170619.1681.100_A', '777962.170619.213.129_A', '777962.170619.241.102_A', '777962.170619.26.103_A', '777962.170619.2779.200_A', '777962.170619.315.112_A', '777962.170619.40.129_B', '777962.170619.408.550_B', '777962.170619.41.101_B', '777962.170619.411.200_B', '777962.170619.417.100_B', '777962.170619.477.104_A', '777962.170619.552.101_A', '777962.170619.99.103_B']

rule filter_target:
	input:
	#jag = expand("filter_logs/{sample}.txt", sample = REF)
		targ = expand("filter_logs/{sample}.txt", sample = SNP50)
	#targ = expand("correct_sex/{sample}.bed", sample = ['58336.170626.1339.200.A', '58336.170626.16709.200.X', '58336.170626.4359.200.Y', '58336.170626.13134.200.Z'])

#58336.170626.94.200.B removed from target because issues were coming up without any sex chromosomes

map_dict = {'777962':"/CIFS/MUG01_N/taylorjerr/PLINK_FILES/9913_HD_161214.map", '227234':"/CIFS/MUG01_N/taylorjerr/PLINK_FILES/9913_GGPF250_161214.map", '58336':"/CIFS/MUG01_N/taylorjerr/PLINK_FILES/9913_SNP50_161214.map", '139977':"/CIFS/MUG01_N/taylorjerr/PLINK_FILES/9913_GGPHDv3_161214.map", '26504':"/CIFS/MUG01_N/taylorjerr/PLINK_FILES/9913_GGPLDv3_161214.map", '30105':"/CIFS/MUG01_N/taylorjerr/PLINK_FILES/9913_GGPLDV4_161214.map", '76999':"/CIFS/MUG01_N/taylorjerr/PLINK_FILES/9913_GGP90KT_161214.map"}
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
	log:
		"logs/duplicate_filter/{sample}.log"
	output:
		dup = "duplicates_filtered/dup_ids/{sample}.txt",
		ped = "duplicates_filtered/{sample}.ped",
	shell:
		"python bin/duplicate_filter.py {input.id} {output.dup} {input.ped} {output.ped}"#; python bin/duplicate_logging.py {output.dup} {output.csv}"



#Dictionaries: This rule takes .map and .ped inputs from dictionaries based on specified output 'wildcard' which is the name of the assay.
#PLINK command (freq) -- Determines call rates for each variant in dataset and returns how many times each snp appears in every assay
#Python Script (snp_call_rate viusalization.py) and (snp_removal_prediction.py).  Takes argument and identifies how many snps will be removed given a specific filter call rate.  Visualization script creates histograms with distributions of call rates and exports to png file in snp_stats/ directory
#Input File Location: raw_genotypes/
#Output File Location: snp_stats/
#Output File Types: (.frq, .log, .nosex)
rule variant_stats:
	input:
		map = mapdicter,
		dup = "duplicates_filtered/dup_ids/{sample}.txt",
		input = "duplicates_filtered/{sample}.ped",
	threads : 4
	params:
		inprefix = "duplicates_filtered/{sample}",
		oprefix = "snp_stats/{sample}"
	benchmark:
		"filter_benchmarks/variant_stats/{sample}.txt"
	output:
		frq = "snp_stats/{sample}.frq",
		log = "snp_stats/{sample}.log",
	shell:
		"plink --file {params.inprefix} --map {input.map} --keep-allele-order --cow --nonfounders --freq --out {params.oprefix}"

#Python Script for Visualization (snp_call_rate_visualization.py) reads in .frq file generated by PLINK's --freq function
#Creates a histogram of the missing genotype rate for each SNP in the dataset.  This call rate is calculated for each SNP as (NCHROBS)/max(NCHROBS)
#Writes figure to a .png file
#Histogram uses 100 bins, and a y-axis max of 5000 SNPs, which look good for existing assays, but can be easily adjusted.
#Input File Location: snp_stats/
#Output File Location: snp_stats/figures/
#Output File Type: (.png)

rule snp_call_rate_visualization:
	input:
		frq = "snp_stats/{sample}.frq"
	output:
		png = "snp_stats/figures/{sample}.snp_call_rate.png"
	shell:
		"python bin/snp_call_rate_visualization.py {input.frq} {output.png}"


#Dictionaries for selecting correct .map and .ped files:  This uses two functions, mapdicter and sampdicter to match the correct .ped and .map files with each assay. Using a specific 'wildcard' for the output of this rule will call the correct input map and ped file for the filter to be used on.
#PLINK filter (mind) -- Filters individuals based on specified call rate.  Individuals missing more than given proportion of their genotypes will be filtered out of the dataset.
#Python script (snp_filtered_log_parsing.py) -- parses output logfile for important metadata and saves to assay's specific csv in /filter_logs directory.
#Input File Location: /raw_genotypes/
#Output File Location: /snp_filtered
#Output: output file suffixes will be (.bed, .bim, .fam, .irem, .log, .nosex)

rule filter_variants:
	input:
		map = mapdicter,
		dup = "duplicates_filtered/dup_ids/{sample}.txt",
		dupped = "duplicates_filtered/{sample}.ped",
		stats = "snp_stats/{sample}.frq",
		png = "snp_stats/figures/{sample}.snp_call_rate.png"
	threads : 4
	params:
		inprefix = "duplicates_filtered/{sample}",
		oprefix="snp_filtered/{sample}",
		logprefix="filter_logs/{sample}"
	benchmark:
		"filter_benchmarks/filter_variants/{sample}.txt"
	output:
		bed="snp_filtered/{sample}.bed",
		bim="snp_filtered/{sample}.bim",
		fam="snp_filtered/{sample}.fam",
		log="snp_filtered/{sample}.log"
	shell:
		"plink --file {params.inprefix} --map {input.map} --threads 4 --keep-allele-order --cow --not-chr 0 --exclude maps/map_issues/snp_ids_to_exclude.txt --geno .1 --make-bed --out {params.oprefix}"#; python bin/snp_filtered_log_parsing.py {output.log} {params.csv}"

rule convert_chip2seq:
		input:
			bed="snp_filtered/{sample}.bed",
			bim="snp_filtered/{sample}.bim",
			fam="snp_filtered/{sample}.fam",
			refalt = "ref_alt/update_all_alleles_170531.txt",
			log="snp_filtered/{sample}.log"
		params:
				iprefix="snp_filtered/{sample}",
				oprefix = "ref_alt/{sample}"
		threads : 4
		benchmark:
				"benchmarks/convert_chip2seq/{sample}.txt"
		log:
				"logs/convert_chip2seq/{sample}.log"
		output:
				bed="ref_alt/{sample}.bed",
				bim="ref_alt/{sample}.bim",
				fam="ref_alt/{sample}.fam",
				log="ref_alt/{sample}.log"
		shell:
				"(plink --bfile {params.iprefix} --cow  --threads 4 --out {params.oprefix} --exclude ref_alt/no_ref_alt_170603.txt  --update-alleles {input.refalt}  --make-bed)> {log}"

rule update_ref:
		input:
			rules.convert_chip2seq.output,
			updateref = "ref_alt/update_all_refs_170531.txt"
		params:
				iprefix="ref_alt/{sample}",
				oprefix = "ref_set/{sample}"
		benchmark:
				"benchmarks/update_ref/{sample}.txt"
		log:
				"logs/update_ref/{sample}.log"
		output:
				bed="ref_set/{sample}.bed",
				bim="ref_set/{sample}.bim",
				fam="ref_set/{sample}.fam",
				log="ref_set/{sample}.log"
		shell:
				"(plink --bfile {params.iprefix} --cow  --out {params.oprefix} --a1-allele {input.updateref}  --make-bed)> {log}"




#PLINK statistic (missing) -- Gathers statistics on individual call rates (reports proportion of missing genotypes for each individual)
#Input File Location: /snp_filtered/
#Output: output file suffixes will be (.imiss, .lmiss, .log, .nosex)
rule individual_stats:
	input:
		rules.update_ref.output,
		bed="ref_set/{sample}.bed"
	params:
		inprefix="ref_set/{sample}",
		oprefix="individual_stats/{sample}"
	benchmark:
		"filter_benchmarks/individual_stats/{sample}.txt"
	output:
		imiss="individual_stats/{sample}.imiss",
		lmiss="individual_stats/{sample}.lmiss",
		log ="individual_stats/{sample}.log"
	shell:
		"plink --bfile {params.inprefix} --cow --missing --real-ref-alleles --out {params.oprefix}" #python bin/individual_call_rate_visualization.py"

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
#Input File Location: /snp_filtered/
#Output File Location: /individual_filtered
#Output: output file suffixes will be (.bed, .bim, .fam, .irem, .log, .nosex)

rule filter_individuals:
	input:
		bed="ref_set/{sample}.bed",
		imiss="individual_stats/{sample}.imiss",
		png="individual_stats/figures/{sample}.individual_call_rate.png"
	params:
		inprefix="ref_set/{sample}",
		oprefix="individual_filtered/{sample}"
	benchmark:
		"filter_benchmarks/filter_individuals/{sample}.txt"
	output:
		bed="individual_filtered/{sample}.bed",
		bim="individual_filtered/{sample}.bim",
		fam="individual_filtered/{sample}.fam",
		log="individual_filtered/{sample}.log"

	shell:
		"plink --bfile {params.inprefix} --cow --mind .1 --real-ref-alleles  --make-bed --out {params.oprefix}"#; python bin/individual_filtered_log_parsing.py {output.log} {params.csv}"



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
		hwe="hwe_stats/{sample}.hwe",
		log="hwe_stats/{sample}.log"
	shell:
		"plink --bfile {params.inprefix}  --real-ref-alleles  --cow --nonfounders --hardy --out {params.oprefix}"

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
#		png="hwe_stats/figures/{sample}.hwe_pvalues.png"
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
		"plink --bfile {params.inprefix} --cow --nonfounders  --real-ref-alleles --hwe 1e-20 --make-bed --out {params.oprefix}"# python bin/hwe_log_parsing.py {output.log} {params.csv}"

#PLINK function (impute-sex) looks at sex provided by ped file, and at X chromosome heterozygosity (and y chromosome variants if provided), and determines whether an animal is male, female, or unknown. If sex from ped file is unknown, this will impute the sex if possible, and write that into the new bed file that it produces.
#Python Scripts: (missexed_animals_filter.py) reads the .sexcheck file produced by impute sex function.  If an individual's sex is changed from known M/F to the opposite sex, it's ID will be written to the {sample}.missexed_animals.txt file, and removed in subsequent step.  (missexed_animals_filter_logging.py) will count lines of filter output, and write to the assay's csv log file.
##Input File Location: hwe_filtered/
##Output File Location: sex_impute/
##Output: output file suffixes will be (.bed, .bim, .fam, .log, .sexcheck, .missexed_animals.txt)
rule impute_sex:
	input:
		bed="hwe_filtered/{sample}.bed"
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
		"plink --bfile {params.inprefix} --cow --real-ref-alleles --impute-sex ycount --make-bed --out {params.oprefix}; python bin/missexed_animals_filter.py {output.sexcheck} {output.txt}"# python bin/missexed_animals_filter_logging.py {output.txt} {output.log} {params.csv}" #{params.logprefix}.csv"

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
		"plink --bfile {params.inprefix} --cow --remove {input.txt} --real-ref-alleles --make-bed --out {params.oprefix}"

rule filter_logging:
	input:
		dup = "duplicates_filtered/dup_ids/{sample}.txt",
		snp = "snp_filtered/{sample}.log",
		ind = "individual_filtered/{sample}.log",
		hwe = "hwe_filtered/{sample}.log",
		sex = "sex_impute/{sample}.log",
		mis = "sex_impute/{sample}.missexed_animals.txt",
		cor = "correct_sex/{sample}.bed"
	benchmark:
		"filter_benchmarks/filter_logging/{sample}.txt"
	log:
		"logs/filter_logging/{sample}.log"
	output:
		log = "filter_logs/{sample}.txt"
	shell:
		"(python bin/combined_filter_logging.py {input.dup} {input.snp} {input.ind} {input.hwe} {input.sex} {input.mis} {output.log}) > {log}"

#This list contains the different file names that will be merged together.
DATA =['227234.170519.1970.GGPF250', '777962.170519.1970.HD']

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
		"python bin/file_list_maker.py correct_sex/allfiles.txt; plink --merge-list correct_sex/allfiles.txt  --cow --real-ref-alleles --make-bed --out {params.oprefix}"



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
		oprefix = "chrsplit/170112_merged.chr{chr}",
		chr = "{chr}"
	benchmark:
		"filter_benchmarks/split_chromosomes/170112_merged.txt"
	log:
		"snake_logs/chromosome_split/170112_merged.chr{chr}.log"
	output:
		bed = "chrsplit/170112_merged.chr{chr}.bed",
		bim = "chrsplit/170112_merged.chr{chr}.bim",
		fam = "chrsplit/170112_merged.chr{chr}.fam",
		log = "chrsplit/170112_merged.chr{chr}.log"
	shell:
		"(plink --bfile {params.inprefix}  --real-ref-alleles --chr {params.chr} --make-bed  --nonfounders --cow --out {params.oprefix})> {log}"



rule assay_chromosomes:
        input:
                bed = expand( "correct_sex/{assay}.bed", assay = DATA),
                bim = expand("correct_sex/{assay}.bim", assay = DATA),
                fam = expand("correct_sex/{assay}.fam", assay = DATA),
                log = expand("correct_sex/{assay}.log", assay = DATA)
        params:
                inprefix = "correct_sex/{sample}",
                oprefix = "eagle_chrsplit/{sample}.chr{chr}",
                chr = "{chr}"
        benchmark:
                "filter_benchmarks/eagle_chrsplit/{sample}.chr{chr}.txt"
        log:
                "snake_logs/eagle_split_chromosomes/{sample}.chr{chr}.log"
        output:
                bed = "assay_chrsplit/{sample}.chr{chr}.bed",
                bim = "assay_chrsplit/{sample}.chr{chr}.bim",
                fam = "assay_chrsplit/{sample}.chr{chr}.fam",
                log = "assay_chrsplit/{sample}.chr{chr}.log"
        shell:
                "(plink --bfile {params.inprefix}  --real-ref-alleles --chr {params.chr} --make-bed  --nonfounders --cow --out {params.oprefix})> {log}"
