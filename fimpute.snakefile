#CHROMOSOMES = list(range(1,31))
DATA =['hol_testset.HD.197','hol_testset.F250.197' , 'hol_testset.SNP50.788','hol_testset.GGPLD.788',]

rule targ:
	input:
		#lege =expand("vcf_to_hap/{sample}.1.chr{chr}.legend.gz",sample = DATA,  chr = list(range(1,30))) #converts phased outputs to haps legend format




#snakemake -s fimpute.snakefile   --cores 8 -np fimpute_per_chip/hol_testset.F250.197.run3.chr25.fimpute
# &> snakerun108.txt
colnums ={'hol_testset.F250.197':2, 'hol_testset.GGPLD.788': 4, 'hol_testset.HD.197' : 1, 'hol_testset.SNP50.788':3}

def colnummer(wildcards):
	return colnums[wildcards.sample]

# need to make an fimpute type file, so any out of step1 needs to go to vcf. This script should work. 
rule vcf_to_fimpute: 
		input:  
			vcf="unphased_vcf/{sample}.run{run}.chr{chr}.vcf",
		params: 
			colnum=colnummer
		log:
			"logs/vcf2fimpute/vcf2impute.run{run}.{chr}.log"
		benchmark:
			"benchmarks/vcf2fimpute/{sample}.run{run}.{chr}.benchmark.txt"
		output: 
			outname="fimpute_per_chip/{sample}.run{run}.chr{chr}.fimpute"
		shell:  
			"python ./bin/vcftoFimpute.py {input.vcf} fimpute {params.colnum} fimpute_per_chip/"

rule fimpute_catter:
	input:
		expand("fimpute_per_chip/{assay}.run{{run}}.chr{{chr}}.fimpute", assay = DATA)
	log:
		"logs/fimpute_catter/{chr}.log"
	benchmark:
		"benchmarks/vcf2fimpute/.run{run}.{chr}.benchmark.txt"
	output :
		"fimpute_genotypes/hol_testset.run{run}.chr{chr}.fimpute"
	shell:
		"(cat fimputehead.txt {input} > {output}) > {log}"
rule fimpute_map_make:
	input:
		expand("assay_chrsplit/{assay}.chr{{chr}}.bim", assay = DATA)
	params:
		chrom = "{chr}"
	log:
		"logs/fimpute_map_make/{chr}.log"
	benchmark:
		"benchmarks/fimpute_map_make/.run{run}.{chr}.benchmark.txt"
	output :
		"fimpute_maps/run{run}.chr{chr}.txt"
	shell:
		"(python ./bin/fimpute_map_maker.py {input} {params.chrom} {output}) > {log}"
		#python ./bin/fimpute_map_maker.py merged_chr 23 fimpute_maps/run4.chr24.txt
		
rule job_writer:
	input:
		"fimpute_ctrs/fimpute_ctr_template.run{run}.txt"
	params:
		chrom = "{chr}"
	output:
		outname = "fimpute_ctrs/{sample}.run{run}.chr{chr}.txt"
	shell:
		"python ./bin/fimpute_job_writer.py {input} {params.chrom} {output.outname}"

rule run_fimpute:
	input:
		gens = "fimpute_genotypes/hol_testset.run{run}.chr{chr}.fimpute",
		ctr = "fimpute_ctrs/{sample}.run{run}.chr{chr}.txt",
		map = "fimpute_maps/run{run}.chr{chr}.txt"
	log:
		"logs/run_fimpute/{run}.{chr}.log"
	benchmark:
		"benchmarks/run_fimpute/run{run}.{chr}.benchmark.txt"
	output :
		"fimpute_output/run{run}.Chr{chr}.{sample}/genotypes_imp.txt"
	shell:
		"(FImpute {input.ctr}) > {log}"


#snakemake -s fimpute.snakefile   --cores 8  fimpute_output/run3.Chr23.hol_testset/genotypes_imp.txt -r  -R fimpute_catter -np



#perl ./bin/vcf2impute_legend_haps.pl -vcf vcf_to_assays/hol_testset.GGPLD.788.1.chr25.vcf -leghap vcf_to_haps/hol_testset.GGPLD.788.1.chr25 -chr 25 



