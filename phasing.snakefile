#CHROMOSOMES = list(range(1,31))
#DATA =['hol_testset.F250.197', 'hol_testset.GGPLD.788', 'hol_testset.HD.197', 'hol_testset.SNP50.788']
DATA =['f250', 'ggpld', 'hd', 'snp50'] #new file names -- these are files that have had ref-alt conversions

rule targ:
	input:
		#jag = expand("eagle_phased_assays/{sample}.chr{chr}.phased.haps.gz",sample = DATA,  chr = list(range(28,31)))
		#jag = expand("merged_chrsplit/hol_testset.merge.1.chr{chr}.bed", chr=list(range(1,31)))
		#jag = expand("vcf_to_assays/{sample}.1.chr{chr}.vcf",sample = DATA,  chr = list(range(22,23))) #phases in parallel then extracts animals on a per assay basis. Ran on 5/23, generated all results
		#lege =expand("vcf_to_hap/{sample}.1.chr{chr}.legend.gz",sample = DATA,  chr = list(range(1,30))) #converts phased outputs to haps legend format
		#jag = "merged_chrsplit/hol_testset.merge.run1.chr1.bed"
		#hapssam =expand("vcf_to_haps/{sample}.run1.chr{chr}.hap.gz",sample = DATA,  chr = list(range(1,30))) #converts phased outputs to haps legend format
		#jag = expand("assay_chrsplit/{assay}.list1.chr{chr}.bed", assay= DATA, chr = list(range(1,30)))
		#jag = expand("eaglemerged/hol_testset.merge.run1.chr{chr}.sample", chr = list(range(1,30)))
		# eaglemerged/hol_testset.merge.run{run}.chr{chr}.sample
		# eaglemerged/hol_testset.merge.run{run}.chr{chr}.sample
		jag = expand("eagle_phased_assays/{sample}.run{run}.chr{chr}.phased.haps.gz", sample = DATA, run = 2, chr = list(range(1,30)))
#snakemake -s phasing.snakefile   --cores 64 -np  &> snakerun108.txt
#snakemake -s phasing.snakefile -p --forceall --cores 64 &> snake_logs/phasing_0603_01.txt ; tail -n 20  snake_logs/phasing_0603_01.txt
hd_or_f250  = {'snp50':"correct_sex/777962.170519.1970.HD",'f250':"correct_sex/227234.170519.1970.GGPF250", 'ggpld':"correct_sex/777962.170519.1970.HD", 'hd':"correct_sex/777962.170519.1970.HD"}
def bedchoice(WC):
	return hd_or_f250[WC.assay]

pull_or_not = {'snp50':"--extract dataprepper/snp50_snps.txt",'f250':"", 'ggpld':"--extract dataprepper/ggpld_snps.txt", 'hd':""}
def snpset(WC):
	return pull_or_not[WC.assay]

run_dict = {'1':'merged_chrsplit', '2':'assay_chrsplit/'}
def runchoice(WC):
	r = WC.run
	samp = WC.sample
	chrom = WC.chr
	if r == '1':
		location = run_dict[r] + '/hol_testset.merge.' + run + '.chr' + chrom +'.bed'
	else:
		location = run_dict[r] + samp + '.list1.chr' + chrom + '.bed'
	return location


rule downsample:
	params:
		idslist = "--keep dataprepper/{assay}_ids.list1.txt",
		extract = snpset,
		bfile = bedchoice,
		oprefix = "testset_assays/{assay}.list1"
	benchmark:
		"benchmarks/downsample/{assay}.list1.txt"
	log:
		"logs/downsample/{assay}.list1.log"
	output:
		bed = "testset_assays/{assay}.list1.bed"
	shell:
		"(plink --bfile {params.bfile}  --real-ref-alleles {params.idslist} {params.extract}  --make-bed  --cow --out {params.oprefix})> {log}"


rule split_chromosomes:
	input:
		bed = expand("testset_assays/{assay}.list1.bed", assay = DATA),
	params:
		inprefix = "testset_assays/{sample}.list1",
		oprefix = "assay_chrsplit/{sample}.list1.chr{chr}",
		chr = "{chr}"
	benchmark:
		"filter_benchmarks/assay_chrsplit/{sample}.chr{chr}.txt"
	log:
		"snake_logs/eagle_split_chromosomes/{sample}.chr{chr}.log"
	output:
		bed = "assay_chrsplit/{sample}.list1.chr{chr}.bed",
		bim = "assay_chrsplit/{sample}.list1.chr{chr}.bim",
		fam = "assay_chrsplit/{sample}.list1.chr{chr}.fam",
		log = "assay_chrsplit/{sample}.list1.chr{chr}.log"
	shell:
		"(plink --bfile {params.inprefix}  --real-ref-alleles --chr {params.chr} --make-bed  --nonfounders --cow --out {params.oprefix})> {log}"

rule run_eagle_single_chrom:
	input:
		bed = runchoice
	params:
		inprefix = "assay_chrsplit/{sample}.list1.chr{chr}",
		oprefix = "eagle_phased_assays/{sample}.run{run}.chr{chr}.phased"
	benchmark:
		"filter_benchmarks/eagle_phased_assays/{sample}.chr{chr}.benchmark.txt"
	threads: 8
	log:
		"snake_logs/eagle_phased_assays/{sample}.chr{chr}.log"
	output:
		sample = "eagle_phased_assays/{sample}.run{run}.chr{chr}.phased.sample",
		haps = "eagle_phased_assays/{sample}.run{run}.chr{chr}.phased.haps.gz"
	shell:
		"(eagle --bfile={params.inprefix} --geneticMapFile=USE_BIM --maxMissingPerSnp 1 --maxMissingPerIndiv 1 --numThreads 8 --outPrefix {params.oprefix}) > {log}"

rule decompress_single_chrom:
		input:
			gzhaps = "eagle_phased_assays/{sample}.run{run}.chr{chr}.phased.haps.gz"
		log:
			"logs/decompress/{sample}.run{run}.chr{chr}.log"
		benchmark:
			"benchmarks/decompress/{sample}.run{run}.chr{chr}.benchmark.txt"
		output:
			haps = "eagle_phased_assays/{sample}.run{run}.chr{chr}.phased.haps"
		shell:
			"(gunzip -c {input.gzhaps} > {output.haps}) > {log}"

rule make_merge_list: #makes a merge list of raw genotypes
	input:
		filelist = expand("assay_chrsplit/{assay}.list1.chr{{chr}}.bed", assay= DATA )
	output:
		"assay_chrsplit/hol_testset.list1.run{run}.chr{chr}.txt"
	shell:
		"python ./bin/merge_file_maker.py {input.filelist} {output}"

rule merged_split_chrs:
	input:
		"assay_chrsplit/hol_testset.list1.run{run}.chr{chr}.txt"
	params:
		outprefix="merged_chrsplit/hol_testset.merge.run{run}.chr{chr}"
	output:
		bed="merged_chrsplit/hol_testset.merge.run{run}.chr{chr}.bed",
		animalset = "merged_chrsplit/hol_testset.merge.run{run}.chr{chr}.fam"
	benchmark:"filter_benchmarks/merged_chrsplit/hol_testset.merg.run{run}.chr{chr}.txt"
	log:
		"logs/merged_chrsplit/hol_testset.merge.run{run}.chr{chr}.log"
	shell:
		"(plink --merge-list {input} --nonfounders --real-ref-alleles --cow --make-bed --out {params.outprefix}) > {log}"
#snakemake -s phasing.snakefile merged_chrsplit/hol_testset.list1.chr22.bed -np --cores 8

#in eagle
#29 1679397 32595 T C
#29 1679398 82560 C T
#29 49415 589629 C A



rule run_eagle_merged:
	input:
		bed = runchoice
		#bed="merged_chrsplit/hol_testset.merge.run{run}.chr{chr}.bed",
	params:
		bed="merged_chrsplit/hol_testset.merge.run{run}.chr{chr}",
		out="eaglemerged/hol_testset.merge.run{run}.chr{chr}"
	threads: 16
	benchmark:
		"benchmarks/eaglemerged/hol_testset.merge.run{run}.chr{chr}.benchmark.txt"
	log:
		"logs/eaglemerged/hol_testset.merge.run{run}.chr{chr}.log"
	output:
		sample = "eaglemerged/hol_testset.merge.run{run}.chr{chr}.sample",
		haps = "eaglemerged/hol_testset.merge.run{run}.chr{chr}.haps.gz"
	shell:
		"(eagle --bfile={params.bed}  --geneticMapFile=USE_BIM --maxMissingPerSnp .99  --maxMissingPerIndiv .99 --numThreads 16 --outPrefix {params.out})> {log} "

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
		haps=temp("eaglemerged/hol_testset.merge.run{run}.chr{chr}.haps"),
	output:
		vcf="makevcf/hol_testset.merge.run{run}.chr{chr}.phased.vcf",
		log="makevcf/logs/hol_testset.merge.run{run}.chr{chr}.log"
	log:
			"logs/makevcf/logs/hol_testset.merge.run{run}.{chr}.log"
	benchmark:
			"benchmarks/makevcf/hol_testset.merge.run{run}.{chr}.benchmark.txt"
	params:
		haps="eaglemerged/hol_testset.merge.run{run}.chr{chr}"
	shell:
		"(shapeit -convert --input-haps {params.haps} --output-log {output.log} --output-vcf {output.vcf}) > {log}"

rule bgzip_vcf:
	input:
		vcf="makevcf/hol_testset.merge.run{run}.chr{chr}.phased.vcf",
	output:
		vcfgz="makevcf/hol_testset.merge.run{run}.chr{chr}.phased.vcf.gz",
		index="makevcf/hol_testset.merge.run{run}.chr{chr}.phased.vcf.gz.tbi",
	log:
		"logs/bgzipvcf/hol_testset.merge.run{run}.chr{chr}"
	benchmark:
		"benchmarks/bgzipvcf/hol_testset.merge.run{run}.chr{chr}"
	shell:
		"(bgzip {input.vcf}; tabix {output.vcfgz}) > {log};"

rule make_vcf_extract_lists:
	input:
		animalset = "merged_chrsplit/hol_testset.merge.run{run}.chr25.fam"
	output:
		keep_ids = "merged_chrsplit/phased_{assay}.list1.run{run}.keepvcf",
		keep_maps = "merged_chrsplit/phased_{assay}.run{run}.vcfregion",
	benchmark:
		"benchmarks/make_vcf_extract_lists/{assay}.run{run}.benchmark.txt"
	log:
		"logs/make_vcf_extract_lists/{assay}.run{run}.benchmark.txt"
	shell:
		"python ./bin/vcfextraction_for_joint_phase.py {input.animalset}"


rule vcf_to_assays: #filter the vcfs on a per assay basis
	input:
		vcfgz="makevcf/hol_testset.merge.run{run}.chr{chr}.phased.vcf.gz",
		index="makevcf/hol_testset.merge.run{run}.chr{chr}.phased.vcf.gz.tbi",
		keep_ids = "merged_chrsplit/phased_{assay}.list1.run{run}.keepvcf",
		keep_maps = "merged_chrsplit/phased_{assay}.run{run}.vcfregion"
	benchmark:
		"benchmarks/vcf_to_assays/{assay}.run{run}.chr{chr}.benchmark.txt"
	log:
		"logs/vcf_to_assays/{assay}.run{run}.chr{chr}.log"
	output:
		vcf = "vcf_to_assays/{assay}.run{run}.chr{chr}.vcf",
	shell:
		"(bcftools view {input.vcfgz} -R {input.keep_maps}  -S {input.keep_ids} -o {output.vcf}) > {log}"

rule vcf_to_hap:
	input:
		vcf = "vcf_to_assays/{assay}.run{run}.chr{chr}.vcf",
	params:
		chr = "{chr}",
		oprefix ="vcf_to_hap/{assay}.run{run}.chr{chr}"
	benchmark:
		"benchmarks/vcf_to_hap/{assay}.run{run}.chr{chr}.benchmark.txt"
	log:
		"logs/vcf_to_haps/{assay}.run{run}.chr{chr}.log"
	output:
		legend = "vcf_to_hap/{assay}.run{run}.chr{chr}.phased.legend",
		hap = "vcf_to_hap/{assay}.run{run}.chr{chr}.phased.haplotypes",
		sample = "vcf_to_hap/{assay}.run{run}.chr{chr}.phased.samples"
	shell:
		"(bcftools convert {input.vcf} --haplegendsample {output.hap}, {output.legend}, {output.sample}) > {log}" #Updated these to use the naming conventions that the shapeit tool outputs the run2 hap,leg,samples, so naming is consistent in Impute2.

rule vcf_to_haps: #doesn't approrpriately name "haps" haps
	input:
		vcf = "vcf_to_assays/{assay}.run{run}.chr{chr}.vcf",
	params:
		chr = "{chr}",
		oprefix ="vcf_to_haps/{assay}.run{run}.chr{chr}.phased"
	benchmark:
		"benchmarks/vcf_to_haps/{assay}.run{run}.chr{chr}.benchmark.txt"
	log:
		"logs/vcf_to_haps/{assay}.run{run}.chr{chr}.log"
	output:
		hap = "vcf_to_haps/{assay}.run{run}.chr{chr}.phased.haps",
		sample = "vcf_to_haps/{assay}.run{run}.chr{chr}.phased.sample"
	shell:
		"(bcftools convert   {input.vcf} --hapsample {output.hap}, {output.sample} ) > {log}"
		#bcftools convert vcf_to_assays/hol_testset.F250.197.1.chr22.vcf --hapsample vcf_to_haps/test.del #produces .hap and .sample
		#bcftools convert  --hapsample vcf_to_haps/test.del
		#"(perl ./bin/vcf2impute_legend_haps.pl -vcf {input.vcf} -leghap {params.oprefix} -chr {params.chr}) 0 >{log}"


#perl ./bin/vcf2impute_legend_haps.pl -vcf vcf_to_assays/hol_testset.GGPLD.788.1.chr25.vcf -leghap vcf_to_haps/hol_testset.GGPLD.788.1.chr25 -chr 25

rule hap_leg:
	input:
		haps = "eagle_phased_assays/{sample}.run{run}.chr{chr}.phased.haps",
		sample = "eagle_phased_assays/{sample}.run{run}.chr{chr}.phased.sample"
	params:
		inprefix = "eagle_phased_assays/{sample}.run{run}.chr{chr}.phased",
		oprefix = "impute_input/{sample}.run{run}.chr{chr}.phased"
	log:
		"snake_logs/hap_leg/{sample}.run{run}.chr{chr}.phased.log"
	benchmark:
		"filter_benchmarks/hap_leg/{sample}.run{run}.chr{chr}.phased.benchmark.txt"
	output:
		hap = "impute_input/{sample}.run{run}.chr{chr}.phased.haplotypes",
		leg = "impute_input/{sample}.run{run}.chr{chr}.phased.legend",
		log = "impute_input/logs/{sample}.run{run}.chr{chr}.phased.log"
	shell:
		"(shapeit -convert --input-haps {params.inprefix} --output-log {output.log} --output-ref {params.oprefix}) > {log}"
