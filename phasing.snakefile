#CHROMOSOMES = list(range(1,31))
DATA =['hol_testset.F250.197', 'hol_testset.GGPLD.788', 'hol_testset.HD.197', 'hol_testset.SNP50.788']

rule targ:
	input:
		#jag = expand("eagle_phased_assays/{sample}.chr{chr}.phased.haps.gz",sample = DATA,  chr = list(range(28,31)))
		#jag = expand("merged_chrsplit/hol_testset.merge.1.chr{chr}.bed", chr=list(range(1,31)))
		#jag = expand("vcf_to_assays/{sample}.1.chr{chr}.vcf",sample = DATA,  chr = list(range(22,23))) #phases in parallel then extracts animals on a per assay basis. Ran on 5/23, generated all results
		#lege =expand("vcf_to_hap/{sample}.1.chr{chr}.legend.gz",sample = DATA,  chr = list(range(1,30))) #converts phased outputs to haps legend format
		hapssam =expand("vcf_to_haps/{sample}.1.chr{chr}.hap.gz",sample = DATA,  chr = list(range(29,30))) #converts phased outputs to haps legend format

rule split_chromosomes:
	input:
		bed = expand( "dataprepper/testset_data/{assay}.bed", assay = DATA),
		bim = expand("dataprepper/testset_data/{assay}.bim", assay = DATA),
		fam = expand("dataprepper/testset_data/{assay}.fam", assay = DATA),
		log = expand("dataprepper/testset_data/{assay}.log", assay = DATA)
	params:
		inprefix = "dataprepper/testset_data/{sample}",
		oprefix = "assay_chrsplit/{sample}.chr{chr}",
		chr = "{chr}"
	benchmark:
		"filter_benchmarks/assay_chrsplit/{sample}.chr{chr}.txt"
	log:
		"snake_logs/eagle_split_chromosomes/{sample}.chr{chr}.log"
	output:
		bed = "assay_chrsplit/{sample}.chr{chr}.bed",
		bim = "assay_chrsplit/{sample}.chr{chr}.bim",
		fam = "assay_chrsplit/{sample}.chr{chr}.fam",
		log = "assay_chrsplit/{sample}.chr{chr}.log"
	shell:
		"(plink --bfile {params.inprefix}  --keep-allele-order --chr {params.chr} --make-bed  --nonfounders --cow --out {params.oprefix})> {log}"

rule run_eagle_single_chrom:
	input:
		bed = "assay_chrsplit/{sample}.chr{chr}.bed",
		bim = "assay_chrsplit/{sample}.chr{chr}.bim",
		fam = "assay_chrsplit/{sample}.chr{chr}.fam",
		log = "assay_chrsplit/{sample}.chr{chr}.log"
	params:
		inprefix = "assay_chrsplit/{sample}.chr{chr}",
		oprefix = "eagle_phased_assays/{sample}.chr{chr}.phased"
	benchmark:
		"filter_benchmarks/eagle_phased_assays/{sample}.chr{chr}.benchmark.txt"
	threads: 8
	log:
		"snake_logs/eagle_phased_assays/{sample}.chr{chr}.log"
	output:
		sample = "eagle_phased_assays/{sample}.chr{chr}.phased.sample",
		haps = "eagle_phased_assays/{sample}.chr{chr}.phased.haps.gz"
	shell:
		"(eagle --bfile={params.inprefix} --geneticMapFile=USE_BIM --maxMissingPerSnp 1 --maxMissingPerIndiv 1 --numThreads 8 --outPrefix {params.oprefix}) > {log}"




rule eagle_to_vcf:
	input:
		haps = "eagle_phased/170112_merged.chr{chr}.haps",
		sample = "eagle_phased/170112_merged.chr{chr}.sample"
	params:
		inprefix = "eagle_phased/170112_merged.chr{chr}",
		oprefix = "eagle_vcfs/170112_merged.chr{chr}"
	benchmark:
		"filter_benchmarks/eagle_to_vcf/170112_merged.chr{chr}.benchmark.txt"
	log:
		"snake_logs/eagle_to_vcf/170112_merged.chr{chr}.log"
	output:
		vcf = "eagle_vcfs/170112_merged.chr{chr}.phased.vcf",
		log = "eagle_vcfs/170112_merged.chr{chr}.log"
	shell:
		"(shapeit -convert --input-haps {params.inprefix} --output-log {output.log} --output-vcf {output.vcf}) > {log}"


rule make_merge_list: #makes a merge list of raw genotypes 
	input:
		f250bed="assay_chrsplit/hol_testset.F250.197.chr{chr}.bed", #from Dataprepper notebook
		snp50bed="assay_chrsplit/hol_testset.SNP50.788.chr{chr}.bed",
		hdbed="assay_chrsplit/hol_testset.HD.197.chr{chr}.bed",
		ggpldbed="assay_chrsplit/hol_testset.GGPLD.788.chr{chr}.bed"
	output:
		"assay_chrsplit/hol_testset.list{run}.chr{chr}.txt"
	shell:
		"python ./bin/merge_file_maker.py {input.f250bed} {input.snp50bed} {input.hdbed} {input.ggpldbed} {output}"

rule merged_split_chrs:		
	input:
		"assay_chrsplit/hol_testset.list{run}.chr{chr}.txt"
	params:
		outprefix="merged_chrsplit/hol_testset.merge.{run}.chr{chr}"
	output:
		bed="merged_chrsplit/hol_testset.merge.{run}.chr{chr}.bed"
	benchmark:"filter_benchmarks/merged_chrsplit/hol_testset.merge{run}.chr{chr}.txt"
	log:
		"logs/merged_chrsplit/hol_testset.merge.{run}.chr{chr}.log"
	shell:
		"(plink --merge-list {input} --nonfounders --cow --make-bed --out {params.outprefix}) > {log}"
#snakemake -s phasing.snakefile merged_chrsplit/hol_testset.list1.chr22.bed -np --cores 8




#rule run_eagle_merged:
#	input:
#		bed="merged_chrsplit/hol_testset.merge.{run}.chr{chr}.bed",
#	params:
#		bed="merged_chrsplit/hol_testset.merge.{run}.chr{chr}",
#		out="eaglemerged/hol_testset.merge.{run}.chr{chr}"
#	threads: 16
#	benchmark:
# 		"benchmarks/eaglemerged/hol_testset.merge.{run}.chr{chr}.benchmark.txt"
# 	log:
# 		"logs/eaglemerged/hol_testset.merge.{run}.chr{chr}.log"
# 	output:
# 		sample = "eaglemerged/hol_testset.merge.{run}.chr{chr}.sample",
# 		haps = "eaglemerged/hol_testset.merge.{run}.chr{chr}.haps.gz"
# 	shell:
# 		"(eagle --bfile={params.bed}  --geneticMapFile=USE_BIM --maxMissingPerSnp .99  --maxMissingPerIndiv .99 --numThreads 16 --outPrefix {params.out})> {log} "

rule decompress:
		input:
			gzhaps = "eaglemerged/{sample}.chr{chr}.haps.gz"
		log:
			"logs/decompress/{sample}.chr{chr}.log"
		benchmark:
			"benchmarks/decompress/{sample}.chr{chr}.benchmark.txt"
		output:
			haps = "eaglemerged/{sample}.chr{chr}.haps"
		shell: 
			"(gunzip -c {input.gzhaps} > {output.haps}) > {log}"

rule make_vcf:
	input:
		haps=temp("eaglemerged/hol_testset.merge.{run}.chr{chr}.haps"),
	output:
		vcf="makevcf/hol_testset.merge.{run}.chr{chr}.phased.vcf",
		log="makevcf/logs/hol_testset.merge.{run}.chr{chr}.log"
	log:
			"logs/makevcf/logs/hol_testset.merge.{run}.{chr}.log"
	benchmark:
			"benchmarks/makevcf/hol_testset.merge.{run}.{chr}.benchmark.txt"
	params:
		haps="eaglemerged/hol_testset.merge.{run}.chr{chr}"
	shell:
		"(shapeit -convert --input-haps {params.haps} --output-log {output.log} --output-vcf {output.vcf}) > {log}"

rule bgzip_vcf:
	input:
		vcf="makevcf/hol_testset.merge.{run}.chr{chr}.phased.vcf",
	output:
		vcfgz="makevcf/hol_testset.merge.{run}.chr{chr}.phased.vcf.gz",
		index="makevcf/hol_testset.merge.{run}.chr{chr}.phased.vcf.gz.tbi",
	log:
		"logs/bgzipvcf/hol_testset.merge.{run}.chr{chr}"
	benchmark:
		"benchmarks/bgzipvcf/hol_testset.merge.{run}.chr{chr}"
	shell:
		"(bgzip {input.vcf}; tabix {output.vcfgz}) > {log};"

rule make_vcf_extract_lists:
	input:
		animalset = "merged_chrsplit/hol_testset.merge.{run}.chr25.fam"
	output:
		keep_maps = "merged_chrsplit/phased_{assay}.{run}.keepvcf",
		keep_ids = "merged_chrsplit/phased_{assay}.{run}.vcfregion", 
	benchmark:
		"benchmarks/make_vcf_extract_lists/{assay}.{run}.benchmark.txt"
	log:
		"logs/make_vcf_extract_lists/{assay}.{run}.benchmark.txt"
	shell:
		"python ./bin/vcfextraction_for_joint_phase.py {input.animalset}" 
		
		



rule vcf_to_assays: #filter the vcfs on a per assay basis
	input:
		vcfgz="makevcf/hol_testset.merge.{run}.chr{chr}.phased.vcf.gz",
		index="makevcf/hol_testset.merge.{run}.chr{chr}.phased.vcf.gz.tbi",
		keep_ids = "merged_chrsplit/phased_{assay}.{run}.keepvcf",
		keep_maps = "merged_chrsplit/phased_{assay}.{run}.vcfregion"
	benchmark:
		"benchmarks/vcf_to_assays/{assay}.{run}.chr{chr}.benchmark.txt"
	log:
		"logs/vcf_to_assays/{assay}.{run}.chr{chr}.log"
	output:
		vcf = "vcf_to_assays/{assay}.{run}.chr{chr}.vcf",
	shell:
		"(bcftools view {input.vcfgz} -R {input.keep_maps}  -S {input.keep_ids} -o {output.vcf}) > {log}"

rule vcf_to_hap:
	input:
		vcf = "vcf_to_assays/{assay}.{run}.chr{chr}.vcf",
	params:
		chr = "{chr}",
		oprefix ="vcf_to_hap/{assay}.{run}.chr{chr}"
	benchmark:
		"benchmarks/vcf_to_hap/{assay}.{run}.chr{chr}.benchmark.txt"
	log:
		"logs/vcf_to_haps/{assay}.{run}.chr{chr}.log"
	output:
		legend = "vcf_to_hap/{assay}.{run}.chr{chr}.legend.gz",
		haps = "vcf_to_hap/{assay}.{run}.chr{chr}.hap.gz",
		sample = "vcf_to_hap/{assay}.{run}.chr{chr}.samples"
	shell:
		"(bcftools convert   -h {params.oprefix}   {input.vcf}) > {log}"

rule vcf_to_haps: #doesn't approrpriately name "haps" haps
	input:
		vcf = "vcf_to_assays/{assay}.{run}.chr{chr}.vcf",
	params:
		chr = "{chr}",
		oprefix ="vcf_to_haps/{assay}.{run}.chr{chr}"
	benchmark:
		"benchmarks/vcf_to_haps/{assay}.{run}.chr{chr}.benchmark.txt"
	log:
		"logs/vcf_to_haps/{assay}.{run}.chr{chr}.log"
	output:
		hap = "vcf_to_haps/{assay}.{run}.chr{chr}.hap.gz",
		sample = "vcf_to_haps/{assay}.{run}.chr{chr}.sample"
	shell:
		"(bcftools convert   {input.vcf} --hapsample {params.oprefix} ) > {log}"
		#bcftools convert vcf_to_assays/hol_testset.F250.197.1.chr22.vcf --hapsample vcf_to_haps/test.del #produces .hap and .sample
		#bcftools convert  --hapsample vcf_to_haps/test.del
		#"(perl ./bin/vcf2impute_legend_haps.pl -vcf {input.vcf} -leghap {params.oprefix} -chr {params.chr}) 0 >{log}" 


#perl ./bin/vcf2impute_legend_haps.pl -vcf vcf_to_assays/hol_testset.GGPLD.788.1.chr25.vcf -leghap vcf_to_haps/hol_testset.GGPLD.788.1.chr25 -chr 25 


#snakemake -s phasing.snakefile   --cores 64  &> snakerun108.txt

